(* ================================================================
   vmmc_2d_thesis.wl  —  VMMC for a Condensate Column Geometry
   ================================================================

   Mathematica translation of the Chapter 2 C++ VMMC implementation
   (Chapter2/src/VMMC.cpp + NucleolusModel.cpp).  Designed for:

     1.  Symbolic detailed-balance checking via check.wls.
     2.  Live animation via animate.wls.

   GEOMETRY
   ─────────────────────────────────────────────────────────────────
   State:   flat array of length L² (row-major, 1-indexed).
            State[[s]] = 0 (empty) or k ∈ {1,...,N} (labeled particle).
            L = Round[Sqrt[Length[state]]].

   Site numbering (same as vmmc_2d.wl):
     (r,c) → s = (r−1)*L + c,   r = row (1..L), c = col (1..L).

   Boundary conditions:
     Row / y direction:  PERIODIC (torus).
     Col / x direction:  HARD WALL — col < 1 or col > L is inaccessible.
                         With $tvHardBarrier = True (default), any
                         cluster move that would push a particle outside
                         the column is rejected outright.

   INTERACTIONS
   ─────────────────────────────────────────────────────────────────
   Distance uses minimum-image in row, absolute in col:
     dr = min-image(r2 − r1, L)   dc = c2 − c1
     d² = dr² + dc²

   Coupling ranges (matching C++ weakD1 / Dsq2 / D2 / Dsq5):
     d²=1,2,4,5  →  γ(x_mid) × $pairJ[a,b]

   A SINGLE coupling symbol J<a><b> per type pair is used for all
   interaction distances.  This is compatible with animate.wls
   (which assigns J12, J13, … from the command line) and correctly
   tests the extended-range VMMC link mechanism under the DB checker.

   GRADIENT (optional, off by default)
   ─────────────────────────────────────────────────────────────────
   γ(x_mid) = Clip[x_mid / L, {0,1}]  when $tvGradient = True
            = 1                         when $tvGradient = False

   x_mid = midpoint column of the bond = (c1 + c2)/2.

   When the gradient is active an internal-bond Metropolis filter
   (C++ VMMC.cpp lines 758–817) is applied after cluster growth:
   ΔE_int = E_post − E_pre over all moving-moving pairs;
   accept with min(1, exp(−β ΔE_int)).  For uniform coupling this
   is always a no-op and draws NO random bits.

   CLUSTER GROWTH
   ─────────────────────────────────────────────────────────────────
   Whitelam–Geissler link-weight mechanism (identical to vmmc_2d.wl):

     wFwd = max(1 − exp(β(eInit − eFwd)), 0)
     wRev = max(1 − exp(β(eInit − eRev)), 0)

   Neighbour shell tested for each cluster particle p: union of the
   √5-shells of p, pPost and pRev.  Required because any site within
   √5 of the virtual forward or reverse position can have a changed
   bond energy with p and MUST be link-tested.

   Hard-wall check for new cluster members is performed in Algorithm
   AFTER the cluster is fully built (not inside the cluster builder),
   so no random bits are consumed by a wall rejection.

   Optional features (all off by default for DB checking):
     $tvCutoff   — cluster-size cutoff ⌊1/r⌋  (C++ proposeMove)
     $tvProbSL   — saturated-link probability  (C++ --phi-sl)
     $tvStokes   — Stokes hydrodynamic drag    (C++ accept())

   CHECKER INTERFACE
   ─────────────────────────────────────────────────────────────────
   Identical to vmmc_2d.wl: BitsToState, DisplayState, ValidStateIDs.
   DynamicSymParams returns the J<a><b> symbols needed for each
   component (same format as vmmc_2d.wl).

   References:
     C++ source:  Chapter2/src/VMMC.cpp, NucleolusModel.cpp
     S. Whitelam & P. L. Geissler, JCP 127, 154101 (2007)
   ================================================================ *)


(* ================================================================
   PARAMETERS  (edit these before running)
   ================================================================ *)

(* Hard wall at column edges.
   True  = closed system → DB checkable.
   False = exit-and-replace (simulation only, breaks DB). *)
$tvHardBarrier = True;

(* Linear chemical gradient γ(x) = col/L  (False = uniform, γ=1).
   When True an internal-bond Metropolis filter corrects DB. *)
$tvGradient = False;

(* Cluster-size cutoff ⌊1/r⌋ (r ~ U[0,1]).  Adds BFS depth. *)
$tvCutoff = False;

(* Saturated-link move probability φ_SL  (0 = disabled). *)
$tvProbSL = 0;
$tvSlN0   = 1;    (* type modulus for SL: skip if type mod $tvSlN0 already in cluster *)

(* Stokes hydrodynamic drag  (False = unit diffusion).
   When True draws one RandomReal[] per multi-particle cluster move. *)
$tvStokes     = False;
$tvRefRadius  = 0.5;   (* reference particle radius for Stokes formula *)


(* ================================================================
   SECTION 1 — Bijective integer encoding   (identical to vmmc_2d.wl)
   ================================================================ *)

$cL[L_]        := $cL[L]      = Sum[Binomial[L,k]*k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]   = Sum[$cL[l], {l, 0, L-1}]
$cLNPre[L_,N_] := $cLNPre[L,N]= Sum[Binomial[L,k]*k!, {k, 0, N-1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L-1, r = rank},
    Do[While[Binomial[x,i] > r, x--]; pos[[i]] = x; r -= Binomial[x,i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n-i]; elems = Delete[elems, idx+1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i-1]]; r = Mod[r, Factorial[i-1]];
       AppendTo[perm, elems[[idx+1]]]; elems = Delete[elems, idx+1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L+1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N+1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]]+1]] = perm[[i]], {i, N}]; arr]


(* ================================================================
   SECTION 2 — Column geometry helpers
   ================================================================ *)

$tvRow[s_, L_] := Ceiling[s / L]
$tvCol[s_, L_] := Mod[s-1, L] + 1

(* Apply displacement {dr, dc} from site s.
   Row: periodic.  Col: hard wall — returns None if outside [1,L]. *)
$tvApplyDir[s_, {dr_, dc_}, L_] :=
  With[{newC = $tvCol[s,L] + dc},
    If[newC < 1 || newC > L, None,
       Mod[$tvRow[s,L] + dr - 1, L] * L + newC]]

(* Squared distance: min-image in row, absolute in col. *)
$tvD2[s1_, s2_, L_] :=
  With[{drRaw = Mod[$tvRow[s1,L] - $tvRow[s2,L], L],
         dc   = $tvCol[s1,L] - $tvCol[s2,L]},
    With[{dr = If[2*drRaw > L, drRaw - L, drRaw]},
      dr^2 + dc^2]]

(* All actual grid sites within 1 ≤ d² ≤ 5 of virtual position {vr, vc}.
   vc may lie outside [1,L] (reverse-virtual move near the wall).
   Uses Select so it returns {} safely when nothing qualifies. *)
$tvVirtualNbrs[{vr_, vc_}, L_] :=
  Select[Range[L^2], Function[q,
    With[{qr = Ceiling[q/L], qc = Mod[q-1,L]+1},
      With[{drRaw = Mod[qr - vr, L], dc = qc - vc},
        With[{dr = If[2*drRaw > L, drRaw - L, drRaw]},
          1 <= dr^2 + dc^2 <= 5]]]]]

(* √5-neighbours of site s (memoised per {s,L}). *)
$tvNbrs[s_, L_] := $tvNbrs[s,L] =
  $tvVirtualNbrs[{$tvRow[s,L], $tvCol[s,L]}, L]


(* ================================================================
   SECTION 3 — Coupling constants and pair energy
   ================================================================ *)

(* Distance-specific coupling symbol per type pair.
   d²=1 → Jd1<lo><hi>,  d²=2 → Jdsq2<lo><hi>,
   d²=4 → Jd2<lo><hi>,  d²=5 → Jdsq5<lo><hi>.
   Each distance range has an independent coupling, matching the
   C++ weakD1 / Dsq2 / D2 / Dsq5 interaction parameters.
   animate.wls assigns these via Jd112=-1.0 Jdsq212=-0.5 etc. *)
$pairJD2[d2_, a_, b_] :=
  If[a == 0 || b == 0, 0,
    With[{lo = Min[a,b], hi = Max[a,b],
          prefix = Switch[d2, 1,"Jd1", 2,"Jdsq2", 4,"Jd2", 5,"Jdsq5", _,None]},
      If[prefix === None, 0,
         ToExpression[prefix <> ToString[lo] <> ToString[hi]]]]]

(* DynamicSymParams: called once per connected component by the checker.
   Returns all four coupling families for the types present. *)
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b,
          {$pairJD2[1,a,b], $pairJD2[2,a,b], $pairJD2[4,a,b], $pairJD2[5,a,b]},
          Nothing],
        {a, types}, {b, types}]|>]

(* AnimationSetupCouplings: hook called by animate.wls after it sets any
   explicitly supplied J symbols.  Fills in missing distance-family symbols
   for each type pair, falling back to J<a><b> if supplied, else random. *)
AnimationSetupCouplings[jArgs_Association, types_List] :=
  Module[{fams = {"Jd1", "Jdsq2", "Jd2", "Jdsq5"}},
    Do[If[a < b,
      With[{baseKey = "J" <> ToString[a] <> ToString[b]},
        Scan[Function[fam,
          With[{key = fam <> ToString[a] <> ToString[b]},
            With[{sym = ToExpression[key]},
              If[!NumericQ[sym],
                Set[Evaluate[sym],
                  Lookup[jArgs, key,
                    Lookup[jArgs, baseKey,
                      RandomVariate[NormalDistribution[-1, 0.3]]]]]]]]],
          fams]]],
      {a, types}, {b, types}]]

(* Gradient scale at bond midpoint column xmid. *)
$tvGamma[xmid_, L_] :=
  If[$tvGradient, Clip[xmid / L, {0, 1}], 1]

(* Pair energy: typeV at virtual position {vr,vc} vs typeQ at site q.
   d²=0 → Infinity (hard core).  d²>5 → 0 (beyond range).
   d²∈{1,2,4,5} → γ(x_mid) × $pairJD2[d², typeV, typeQ]. *)
$tvVPairE[typeV_, {vr_, vc_}, typeQ_, q_, L_] :=
  If[typeV == 0 || typeQ == 0, 0,
    With[{qr = Ceiling[q/L], qc = Mod[q-1,L]+1},
      With[{drRaw = Mod[qr - vr, L], dc = qc - vc},
        With[{dr = If[2*drRaw > L, drRaw - L, drRaw],
               g = $tvGamma[(vc + qc)/2, L]},
          With[{d2 = dr^2 + dc^2},
            Which[
              d2 == 0, Infinity,
              d2 > 5,  0,
              True,    g * $pairJD2[d2, typeV, typeQ]]]]]]]


(* ================================================================
   SECTION 4 — Energy function
   ================================================================ *)

(* Unique undirected bonds at 1 ≤ d² ≤ 5 on the L×L column grid.
   Stored as {s1, s2, d²} triples so energy[] can look up the correct
   coupling symbol without recomputing distances.  Memoised per L. *)
$tvUniqueBonds[L_] := $tvUniqueBonds[L] =
  Flatten[
    Table[With[{d2 = $tvD2[s1,s2,L]},
            If[1 <= d2 <= 5, {{s1, s2, d2}}, {}]],
          {s1, L^2}, {s2, s1+1, L^2}], 2]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total @ Map[
      Function[bond,
        With[{a = state[[bond[[1]]]], b = state[[bond[[2]]]]},
          If[a == 0 || b == 0, 0,
            $tvGamma[($tvCol[bond[[1]],L] + $tvCol[bond[[2]],L])/2, L] *
            $pairJD2[bond[[3]], a, b]]]],
      $tvUniqueBonds[L]]]


(* ================================================================
   SECTION 5 — Stokes drag helper  (used when $tvStokes = True)
   ================================================================

   Scale factor = R_ref / (R_ref + sqrt(H/n)) following C++ accept().
   H = sum over cluster of (perpendicular component of delta from CoM)²,
   where perpendicular to dir = {dr,dc} is given by the 2D cross product:
     perp_i = (col_i − colCoM) × dr − (row_i − rowCoM) × dc
   For a single particle H = 0 → scaleFactor = 1 (no rejection).   *)

$tvStokesScale[cluster_, dir_, L_] :=
  With[{n = Length[cluster],
        rows = Map[$tvRow[#, L] &, cluster],
        cols = Map[$tvCol[#, L] &, cluster]},
    With[{rCoM = Mean[rows], cCoM = Mean[cols]},
      With[{perp = (cols - cCoM) * dir[[1]] - (rows - rCoM) * dir[[2]]},
        $tvRefRadius / ($tvRefRadius + Sqrt[Total[perp^2] / n])]]]


(* ================================================================
   SECTION 6 — VMMC cluster builder
   ================================================================

   Structure is IDENTICAL to $vmmcBuildCluster in vmmc_2d.wl; the
   only differences are:
     • √5-neighbour shells (not 1-shells) via $tvVirtualNbrs.
     • Virtual positions stored as {vr,vc} pairs (vc may be off-grid).
     • No hard-wall check inside this function — wall is checked in
       Algorithm after the full cluster is assembled.

   Returns the cluster site list on success, or None on frustrated link.
   ================================================================ *)

$tvBuildCluster[state_, L_, seed_, dir_, slTypes_] :=
  Module[{
    cluster   = {seed},
    inCluster = <|seed -> True|>,
    queue     = {seed},
    frustrated = False,
    p, pType, pRC, pPostRC, pRevRC, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },

    While[queue =!= {} && !frustrated,

      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];

      (* Virtual positions for p: forward (pPost) and reverse (pRev).
         Row wraps periodically; col is unconstrained (may be 0 or L+1).
         These are {vr,vc} pairs, not site indices. *)
      pRC     = {$tvRow[p,L], $tvCol[p,L]};
      pPostRC = {Mod[pRC[[1]] + dir[[1]] - 1, L] + 1, pRC[[2]] + dir[[2]]};
      pRevRC  = {Mod[pRC[[1]] - dir[[1]] - 1, L] + 1, pRC[[2]] - dir[[2]]};

      (* Union of √5-shells of p, pPost, pRev.
         Any site within √5 of any of these three positions may have a
         changed pair energy with p and must receive a link test.
         DeleteDuplicates prevents double-testing (critical for L=2). *)
      nbrs = DeleteDuplicates @ Join[
               $tvNbrs[p, L],
               $tvVirtualNbrs[pPostRC, L],
               $tvVirtualNbrs[pRevRC,  L]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          (* SL gating: skip if this particle's type is already in cluster *)
          If[$tvProbSL > 0 && $tvSlN0 > 1 &&
             KeyExistsQ[slTypes, Mod[qType, $tvSlN0]],
            Continue[]];

          (* Three pair energies for the link-weight formula *)
          eInit = $tvVPairE[pType, pRC,     qType, q, L];
          eFwd  = $tvVPairE[pType, pPostRC, qType, q, L];
          eRev  = $tvVPairE[pType, pRevRC,  qType, q, L];

          (* wFwd = max(1 − exp(β(eInit−eFwd)), 0).
             Piecewise — not Max — so FullSimplify can resolve J sign cases.
             Leading Infinity guards prevent ComplexInfinity with symbolic β.
             (Identical logic to vmmc_2d.wl.) *)
          wFwd = Piecewise[{
              {1,                               eFwd === Infinity},
              {1 - Exp[\[Beta](eInit - eFwd)],  eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                               eRev === Infinity},
              {1 - Exp[\[Beta](eInit - eRev)],  eInit < eRev}},
            0];

          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1,
                      eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta](eInit - eRev)],
                      eFwd === Infinity && eInit < eRev},
                  {0,
                      eFwd === Infinity},
                  {1,
                      eRev === Infinity && eInit < eFwd},
                  {Min[(1-Exp[\[Beta](eInit-eRev)]) /
                        (1-Exp[\[Beta](eInit-eFwd)]), 1],
                      eInit < eFwd && eInit < eRev},
                  {0,
                      eInit < eFwd}},
                0],
              (* Frustrated link: reject entire move *)
              frustrated = True; Break[],

              (* Link accepted: q joins the cluster *)
              AppendTo[cluster, q];
              inCluster[q] = True;
              If[$tvProbSL > 0 && $tvSlN0 > 1,
                AssociateTo[slTypes, Mod[qType, $tvSlN0] -> True]];
              AppendTo[queue, q]
            ]
          ]
        ],
        {k, Length[nbrs]}
      ]
    ];

    If[frustrated, None, cluster]
  ]


(* ================================================================
   SECTION 7 — Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{
    L, occupied, seed, dirs, dir, slTypes, cutOff,
    cluster, newState, dE, dEInt, sA, sB, destA, destB, ePre, ePost
  },

    L = Round[Sqrt[Length[state]]];

    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    (* Four cardinal directions: {dr, dc} = row-change, col-change. *)
    dirs = {{0,1}, {0,-1}, {1,0}, {-1,0}};
    dir  = RandomChoice[dirs];

    (* ── Hard-wall pre-check: seed destination outside column? ────────
       No random bits consumed — pure geometry. *)
    If[$tvHardBarrier && $tvApplyDir[seed, dir, L] === None,
      Return[state]];

    (* ── Optional: cluster-size cutoff  (C++ proposeMove) ────────────
       Draw cutOff = ⌊1/r⌋ before building to match C++ RNG ordering. *)
    cutOff = If[$tvCutoff,
               With[{r = RandomReal[]}, If[r == 0, Infinity, Floor[1/r]]],
               Infinity];

    (* ── Optional: saturated-link mode ───────────────────────────────
       Seed's own type is pre-loaded into slTypes. *)
    slTypes = If[$tvProbSL > 0 && $tvSlN0 > 1 && RandomReal[] < $tvProbSL,
                 <|Mod[state[[seed]], $tvSlN0] -> True|>,
                 <||>];

    (* ── Grow the VMMC cluster ────────────────────────────────────── *)
    cluster = $tvBuildCluster[state, L, seed, dir, slTypes];
    If[cluster === None, Return[state]];    (* frustrated link *)

    (* ── Hard-wall check for all cluster members ──────────────────────
       Done here (after building, before any move) so NO random bits are
       consumed by a wall rejection.  Superdetailed balance is intact
       because this is a deterministic per-cluster-geometry check. *)
    If[$tvHardBarrier &&
       AnyTrue[cluster, $tvApplyDir[#, dir, L] === None &],
      Return[state]];

    (* ── Cluster-size cutoff check ────────────────────────────────── *)
    If[Length[cluster] > cutOff, Return[state]];

    (* ── Optional: Stokes drag rejection ─────────────────────────────
       scaleFactor = R_ref / (R_ref + sqrt(H/n)).  For single particles
       scaleFactor = 1 (never rejected).  Drawn after cutoff for C++ fidelity. *)
    If[$tvStokes && Length[cluster] > 1,
      If[RandomReal[] > $tvStokesScale[cluster, dir, L],
        Return[state]]];

    (* ── Apply rigid cluster translation ─────────────────────────── *)
    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[
      With[{dest = $tvApplyDir[cluster[[i]], dir, L]},
        (* dest = None cannot happen: hard-wall check above ensures safety *)
        If[newState[[dest]] =!= 0, Return[state, Module]];   (* overlap *)
        newState[[dest]] = state[[cluster[[i]]]]],
      {i, Length[cluster]}];

    (* ── Internal-bond Metropolis filter ─────────────────────────────
       Needed for DB when $tvGradient = True:
         ΔE_int = Σ_{a<b in cluster} [E_post(a,b) − E_pre(a,b)]
       The pair energies E_pre/E_post use the γ-scaled coupling.
       For $tvGradient = False:  every ΔE = 0, no random bit consumed.
       For $tvGradient = True:   one RandomReal[] drawn via MetropolisProb. *)
    If[$tvGradient,
      dEInt = Total @ Flatten @ Table[
        Module[{
          sA = cluster[[a]], sB = cluster[[b]],
          rcA, rcB, destA, destB, ePre2, ePost2},
          rcA   = {$tvRow[sA,L], $tvCol[sA,L]};
          rcB   = {$tvRow[sB,L], $tvCol[sB,L]};
          destA = $tvApplyDir[sA, dir, L];
          destB = $tvApplyDir[sB, dir, L];
          ePre2  = $tvVPairE[state[[sA]], rcA,
                             state[[sB]], sB,    L];
          ePost2 = $tvVPairE[state[[sA]], {$tvRow[destA,L],$tvCol[destA,L]},
                             state[[sB]], destB,  L];
          If[ePre2 > 1*^5 || ePost2 > 1*^5, 0, ePost2 - ePre2]
        ],
        {a, Length[cluster]}, {b, a+1, Length[cluster]}];
      If[RandomReal[] >= MetropolisProb[dEInt], Return[state]]];

    (* ── Final energy change (no-op acceptance; VMMC DB is link-enforced)
       MetropolisProb[dE] is called for checker interface compliance but
       no RandomReal[] is drawn here — the result is discarded. *)
    dE = energy[newState] - energy[state];
    MetropolisProb[dE];

    newState
  ]


(* ================================================================
   SECTION 8 — Checker interface   (identical to vmmc_2d.wl)
   ================================================================ *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2+1]-1, maxId]]];
      L++];
    ids]

numBeta = 1
