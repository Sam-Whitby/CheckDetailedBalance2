(* ================================================================
   vmmc_2d_thesis.wl  —  VMMC for the Column / Condensate Geometry
   ================================================================

   Mathematica translation of the Chapter 2 C++ VMMC implementation
   (Chapter2/src/VMMC.cpp + NucleolusModel.cpp).  Designed for:

     1.  Symbolic detailed-balance checking (hard-barrier mode).
     2.  Live animation via animate.wls.

   GEOMETRY
   ─────────────────────────────────────────────────────────────────
   State:   flat array of length L² (row-major, 1-indexed).
            State[[s]] = 0 (empty) or k ∈ {1,...,N} (labeled particle).
            L = Round[Sqrt[Length[state]]].

   Site (r, c):  r = row (y-direction, PERIODIC),  c = col (x-direction).
   Index:        s = (r-1)*L + c,   r,c ∈ {1,...,L}.

   Boundary conditions:
     y / row direction:   periodic (torus)
     x / col direction:   hard wall — col 0 and col L+1 are inaccessible.
                          With $tvHardBarrier = True (default for DB check)
                          any cluster move that would carry a particle
                          outside the column is rejected outright.

   INTERACTIONS
   ─────────────────────────────────────────────────────────────────
   Distance is computed with minimum-image in y and absolute in x:
     dc = c2 − c1  (absolute)
     dr = min-image(r2 − r1, L)   (periodic)
     d² = dr² + dc²

   Four coupling ranges (identical to C++ weakD1/Dsq2/D2/Dsq5):
     d²=1  (cardinal,    distance 1):   γ(x_mid) × Jd1_<a><b>
     d²=2  (diagonal,    distance √2):  γ(x_mid) × Jdsq2_<a><b>
     d²=4  (cardinal×2,  distance 2):   γ(x_mid) × Jd2_<a><b>
     d²=5  (knight move, distance √5):  γ(x_mid) × Jdsq5_<a><b>

   Hard-core exclusion (d²=0, same site) → Infinity.

   GRADIENT
   ─────────────────────────────────────────────────────────────────
   γ(x_mid) = x_mid / L   if  $tvGradient = True
            = 1            if  $tvGradient = False  (default, uniform)
   x_mid = (c1 + c2)/2  (column midpoint of the bond).

   When the gradient is active, the VMMC link weights correctly use
   the position-dependent energy, and an internal-bond Metropolis
   filter (C++ lines 758–817 of VMMC.cpp) is applied after cluster
   growth to restore detailed balance for the co-moving pairs.

   CLUSTER GROWTH
   ─────────────────────────────────────────────────────────────────
   Whitelam–Geissler link-weight mechanism (same as vmmc_2d.wl):

     wFwd = max(1 − exp(β(eInit − eFwd)), 0)
     wRev = max(1 − exp(β(eInit − eRev)), 0)

   Neighbour shell: union of all sites within d=√5 of p, pPost, pRev.
   This is essential for correctness at extended range, as any site
   within √5 of the virtual forward/reverse position can change its
   bond energy and must be subjected to a link test.

   Optional cluster-size cutoff (C++ VMMC.cpp proposeMove):
     draw r ~ U[0,1],  cutOff = ⌊1/r⌋
     reject move if |cluster| > cutOff

   Optional saturated-link (SL) mode (C++ --phi-sl flag):
     with probability $tvProbSL, skip neighbours whose type
     (mod $tvSlN0) is already present in the cluster.

   CHECKER INTERFACE
   ─────────────────────────────────────────────────────────────────
   Compatible with check.wls / dbc_core.wl (same structure as
   vmmc_2d.wl). BitsToState filters to perfect-square arrays.
   Coupling symbols Jd1_ab, Jdsq2_ab, Jd2_ab, Jdsq5_ab replace the
   single J_ab of vmmc_2d.wl.

   References:
     C++ source:  Chapter2/src/VMMC.cpp, NucleolusModel.cpp
     S. Whitelam & P. L. Geissler, JCP 127, 154101 (2007)
   ================================================================ *)


(* ---- Simulation parameters (change before running) ---------------------- *)

(* Hard wall at column boundaries.  True = hard barrier (no exit; system
   closed → detailed balance can be checked).
   False = exit-and-replace condition (breaks DB; for simulation only). *)
$tvHardBarrier = True;

(* Linear chemical gradient γ(x) = col/L.  False = uniform coupling (γ=1).
   When True an internal-bond Metropolis filter is applied after cluster
   growth to ensure detailed balance with spatially-varying energies. *)
$tvGradient = False;

(* Cluster-size cutoff: draw r~U[0,1], reject move if |cluster|>⌊1/r⌋.
   Set False to simplify symbolic DB check (cutoff does not affect DB but
   adds BFS branches). *)
$tvCutoff = False;

(* Saturated-link move probability (0 = disabled).
   When enabled, a fraction $tvProbSL of moves use the SL rule:
   no two cluster members share the same type mod $tvSlN0. *)
$tvProbSL = 0;
$tvSlN0   = 1;


(* ================================================================
   SECTION 1 — Bijective integer encoding   (identical to vmmc_2d.wl)
   ================================================================ *)

$cL[L_]        := $cL[L]     = Sum[Binomial[L, k]*k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]  = Sum[$cL[l], {l, 0, L-1}]
$cLNPre[L_,N_] := $cLNPre[L,N] = Sum[Binomial[L,k]*k!, {k, 0, N-1}]

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

(* Row and column of site s in an L×L grid (1-indexed, row-major).
   Row = y-direction (periodic).  Col = x-direction (hard wall). *)
$tvRow[s_, L_] := Ceiling[s / L]
$tvCol[s_, L_] := Mod[s-1, L] + 1

(* Apply displacement {dr, dc} from site s.
   Row: periodic wrap.  Col: hard wall — returns None if out of [1,L]. *)
$tvApplyDir[s_, {dr_, dc_}, L_] :=
  Module[{newC = $tvCol[s,L] + dc},
    If[newC < 1 || newC > L, None,
       (Mod[$tvRow[s,L] + dr - 1, L]) * L + newC]]

(* Distance squared between s1 and s2: absolute x, minimum-image y. *)
$tvD2[s1_, s2_, L_] :=
  Module[{dr = Mod[$tvRow[s1,L] - $tvRow[s2,L], L],
          dc = $tvCol[s1,L] - $tvCol[s2,L]},
    If[2*dr > L, dr -= L];
    dr^2 + dc^2]

(* All actual grid sites within 1 ≤ d² ≤ 5 of virtual position {vr, vc}.
   vc may be outside [1,L] (e.g. for a reverse-virtual move at the wall).
   Row is periodic; col is bounded to [1,L]. *)
$tvVirtualNbrs[{vr_, vc_}, L_] :=
  Reap[
    Do[
      If[1 <= qc <= L,
        Module[{dr = Mod[qr - vr, L], dc = qc - vc, d2},
          If[dr > L/2, dr -= L];
          d2 = dr^2 + dc^2;
          If[1 <= d2 <= 5, Sow[(qr-1)*L + qc]]]],
      {qr, 1, L}, {qc, 1, L}
    ]][[2, 1]] /. {} -> {}      (* ensure list even when empty *)

(* All sites within d²≤5 of site s (cached per site). *)
$tvNbrs[s_, L_] := $tvNbrs[s,L] =
  $tvVirtualNbrs[{$tvRow[s,L], $tvCol[s,L]}, L]


(* ================================================================
   SECTION 3 — Coupling constants and pair energy
   ================================================================ *)

(* Coupling symbol for type-pair {a,b} at distance label dist.
   Returns 0 if either type is 0 (hole).  Names:
     Jd1_<min><max>   (d²=1)
     Jdsq2_<min><max> (d²=2)
     Jd2_<min><max>   (d²=4)
     Jdsq5_<min><max> (d²=5) *)
$tvSym[a_, b_, prefix_String] :=
  If[a == 0 || b == 0, 0,
    ToExpression[prefix <> ToString[Min[a,b]] <> ToString[Max[a,b]]]]

$tvJd1[a_,b_]   := $tvSym[a, b, "Jd1"]
$tvJdsq2[a_,b_] := $tvSym[a, b, "Jdsq2"]
$tvJd2[a_,b_]   := $tvSym[a, b, "Jd2"]
$tvJdsq5[a_,b_] := $tvSym[a, b, "Jdsq5"]

(* Gradient scale factor at bond-midpoint column xmid. *)
$tvGamma[xmid_, L_] :=
  If[$tvGradient, Clip[xmid / L, {0, 1}], 1]

(* Pair energy between a particle of type typeV at virtual position
   {vr, vc} (vc may be outside [1,L]) and a particle of type typeQ
   at actual site q.  Same-site → Infinity (hard core). *)
$tvVPairE[typeV_, {vr_, vc_}, typeQ_, q_, L_] :=
  Module[{qr = $tvRow[q,L], qc = $tvCol[q,L], dr, dc, d2, g},
    If[typeV == 0 || typeQ == 0, Return[0]];
    dc = qc - vc;
    dr = Mod[qr - vr, L];
    If[2*dr > L, dr -= L];
    d2 = dr^2 + dc^2;
    Which[
      d2 == 0, Infinity,   (* same site: hard core *)
      d2 > 5,  0,          (* beyond √5 range *)
      True,
        g = $tvGamma[(vc + qc)/2, L];
        Switch[d2,
          1, g * $tvJd1[typeV, typeQ],
          2, g * $tvJdsq2[typeV, typeQ],
          4, g * $tvJd2[typeV, typeQ],
          5, g * $tvJdsq5[typeV, typeQ],
          _, 0]]]

(* DynamicSymParams: returns all symbolic coupling parameters needed for
   the given component states (called by the DB checker). *)
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a <= b,
          {$tvSym[a, b, "Jd1"],
           $tvSym[a, b, "Jdsq2"],
           $tvSym[a, b, "Jd2"],
           $tvSym[a, b, "Jdsq5"]},
          Nothing],
        {a, types}, {b, types}]|>]


(* ================================================================
   SECTION 4 — Energy function
   ================================================================ *)

(* All unique bonds on the L×L column grid at 1 ≤ d² ≤ 5.
   Memoised per L. *)
$tvUniqueBonds[L_] := $tvUniqueBonds[L] =
  Flatten[
    Table[
      If[1 <= $tvD2[s1, s2, L] <= 5, {{s1, s2}}, {}],
      {s1, L^2}, {s2, s1+1, L^2}],
    2]

(* Total energy of a state.  All couplings symbolic; positions numeric. *)
energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total @ Map[
      Function[bond,
        Module[{a = state[[bond[[1]]]], b = state[[bond[[2]]]], d2, g},
          If[a == 0 || b == 0, 0,
            d2 = $tvD2[bond[[1]], bond[[2]], L];
            g  = $tvGamma[($tvCol[bond[[1]],L] + $tvCol[bond[[2]],L])/2, L];
            Switch[d2,
              1, g * $tvJd1[a,b],
              2, g * $tvJdsq2[a,b],
              4, g * $tvJd2[a,b],
              5, g * $tvJdsq5[a,b],
              _, 0]]]],
      $tvUniqueBonds[L]]]


(* ================================================================
   SECTION 5 — VMMC cluster builder
   ================================================================

   Identical structure to vmmc_2d.wl's $vmmcBuildCluster, extended to:
     • Union of √5-neighbour shells of p, pPost, pRev.
     • Hard-barrier check: any cluster particle whose forward-virtual
       site is outside the column → reject immediately (return None).
     • Optional cluster-size cutoff (drawn before building; checked after).
     • Optional saturated-link (SL) recruitment gating.
     • Internal-bond Metropolis filter applied in Algorithm[] after
       building (needed when $tvGradient = True).

   Returns a list of cluster site indices on success, or None if any
   frustrated link is encountered or the hard barrier is violated.
   ================================================================ *)

$tvBuildCluster[state_, L_, seed_, dir_, slTypes_] :=
  Module[{
    cluster    = {seed},
    inCluster  = <|seed -> True|>,
    queue      = {seed},
    frustrated = False,
    p, pType, pRC, pPostRC, pRevRC, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },

    (* ── Hard barrier: check seed's forward virtual position first ── *)
    If[$tvHardBarrier && $tvApplyDir[seed, dir, L] === None,
      Return[None]];

    While[queue =!= {} && !frustrated,

      p      = First[queue]; queue = Rest[queue];
      pType  = state[[p]];

      (* Virtual (row, col) for forward and reverse moves.
         Col may go outside [1,L] for reverse moves near the wall;
         that is intentional — pair energies are still computed
         (matching C++ behaviour) while the hard barrier blocks actual
         moves that would go outside. *)
      pRC     = {$tvRow[p,L], $tvCol[p,L]};
      pPostRC = {Mod[pRC[[1]] + dir[[1]] - 1, L] + 1,
                 pRC[[2]] + dir[[2]]};
      pRevRC  = {Mod[pRC[[1]] - dir[[1]] - 1, L] + 1,
                 pRC[[2]] - dir[[2]]};

      (* Union of √5-neighbour shells of p, pPost, pRev.
         Handles L=2 wrap-arounds via DeleteDuplicates. *)
      nbrs = DeleteDuplicates @ Join[
               $tvNbrs[p, L],
               $tvVirtualNbrs[pPostRC, L],
               $tvVirtualNbrs[pRevRC,  L]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          (* SL gating: skip neighbour if its type already in cluster *)
          If[$tvProbSL > 0 && $tvSlN0 > 1 &&
             MemberQ[Keys[slTypes], Mod[qType, $tvSlN0]],
            Continue[]];

          (* Pair energies at current, forward-virtual, reverse-virtual *)
          eInit = $tvVPairE[pType, pRC,     qType, q, L];
          eFwd  = $tvVPairE[pType, pPostRC, qType, q, L];
          eRev  = $tvVPairE[pType, pRevRC,  qType, q, L];

          (* Forward link weight  wFwd = max(1 − exp(β(eInit−eFwd)), 0)
             Piecewise (not Max) so FullSimplify can resolve sign cases.
             Leading Infinity guards prevent ComplexInfinity with symbolic β. *)
          wFwd = Piecewise[{
              {1,                                eFwd === Infinity},
              {1 - Exp[\[Beta] (eInit - eFwd)],  eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                                eRev === Infinity},
              {1 - Exp[\[Beta] (eInit - eRev)],  eInit < eRev}},
            0];

          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1, eFwd===Infinity && eRev===Infinity},
                  {1 - Exp[\[Beta](eInit-eRev)],
                      eFwd===Infinity && eInit < eRev},
                  {0, eFwd===Infinity},
                  {1, eRev===Infinity && eInit < eFwd},
                  {Min[(1-Exp[\[Beta](eInit-eRev)]) /
                       (1-Exp[\[Beta](eInit-eFwd)]), 1],
                      eInit < eFwd && eInit < eRev},
                  {0, eInit < eFwd}},
                0],
              (* Frustrated link → reject entire move *)
              frustrated = True; Break[],

              (* Hard barrier: reject if new member q can't move forward *)
              If[$tvHardBarrier && $tvApplyDir[q, dir, L] === None,
                frustrated = True; Break[]];

              AppendTo[cluster, q];
              inCluster[q] = True;
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
   SECTION 6 — Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{
    L, occupied, seed, dirs, dir, isSL, slTypes,
    cutOff, cluster, newState, dE, dEInt,
    clusterRC, newRC, a, b, ePre, ePost
  },

    L = Round[Sqrt[Length[state]]];

    (* Select a random occupied site as the seed *)
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    (* Choose one of the four cardinal directions.
       Only cardinal moves are used (matching C++ lattice VMMC with 4 neighbours).
       Diagonal / extended-range moves are handled via the link weights. *)
    dirs = {{0,1}, {0,-1}, {1,0}, {-1,0}};    (* right, left, down, up *)
    dir  = RandomChoice[dirs];

    (* ── Optional: Saturated-link mode ───────────────────────────── *)
    isSL = ($tvProbSL > 0) && (RandomReal[] < $tvProbSL);
    slTypes = If[isSL && $tvSlN0 > 1,
                 <|Mod[state[[seed]], $tvSlN0] -> True|>,
                 <||>];

    (* ── Optional: cluster-size cutoff (draw before building) ─────── *)
    cutOff = If[$tvCutoff,
               Module[{r = RandomReal[]},
                 If[r == 0, Infinity, Floor[1/r]]],
               Infinity];

    (* ── Build the VMMC cluster ───────────────────────────────────── *)
    cluster = $tvBuildCluster[state, L, seed, dir, slTypes];
    If[cluster === None, Return[state]];    (* frustrated or boundary *)

    (* ── Apply cluster-size cutoff ────────────────────────────────── *)
    If[Length[cluster] > cutOff, Return[state]];

    (* ── Apply the rigid cluster translation ─────────────────────── *)
    newState = state;
    Do[newState[[cluster[[i]]]] = 0,       {i, Length[cluster]}];
    Do[
      With[{dest = $tvApplyDir[cluster[[i]], dir, L]},
        If[dest === None, Return[state, Module]];         (* should not happen: hard barrier already checked *)
        If[newState[[dest]] =!= 0, Return[state, Module]]; (* overlap *)
        newState[[dest]] = state[[cluster[[i]]]]
      ],
      {i, Length[cluster]}
    ];

    (* ── Internal-bond Metropolis filter ─────────────────────────────
       Required for detailed balance when interactions are spatially
       varying (gradient active).  Computes ΔE over all moving-moving
       pairs and applies min(1, exp(−βΔE)).
       For uniform coupling ($tvGradient=False) every term is zero
       (translation-invariant energy), so this is always a no-op.
       MetropolisProb is called unconditionally for checker compliance.
       ─────────────────────────────────────────────────────────────── *)
    dEInt = Total @ Flatten @ Table[
      Module[{
        sA = cluster[[a]], sB = cluster[[b]],
        rcA = {$tvRow[cluster[[a]], L], $tvCol[cluster[[a]], L]},
        rcB = {$tvRow[cluster[[b]], L], $tvCol[cluster[[b]], L]},
        rcANew, rcBNew},
        rcANew = {Mod[rcA[[1]] + dir[[1]] - 1, L] + 1, rcA[[2]] + dir[[2]]};
        rcBNew = {Mod[rcB[[1]] + dir[[1]] - 1, L] + 1, rcB[[2]] + dir[[2]]};
        ePre  = $tvVPairE[state[[sA]], rcA,    state[[sB]], sB, L];
        ePost = $tvVPairE[state[[sA]], rcANew, state[[sB]],
                  (* site index of new position of B *)
                  (rcBNew[[1]] - 1)*L + rcBNew[[2]], L];
        If[ePre > 1*^5 || ePost > 1*^5, 0, ePost - ePre]
      ],
      {a, Length[cluster]}, {b, a+1, Length[cluster]}];

    (* MetropolisProb[dEInt] = min(1, exp(−β dEInt)).
       If dEInt=0 (no gradient), this equals 1 and the RandomReal[]
       comparison below is always satisfied (no rejection). *)
    If[RandomReal[] >= MetropolisProb[dEInt], Return[state]];

    (* ── Final energy delta (for interface; VMMC DB is link-enforced) *)
    dE = energy[newState] - energy[state];
    MetropolisProb[dE];   (* called for checker interface compliance *)

    newState
  ]


(* ================================================================
   SECTION 7 — Checker interface   (same structure as vmmc_2d.wl)
   ================================================================ *)

(* Accept only IDs whose decoded array has perfect-square length (L×L). *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

(* Display state as a grid: {1,2}|{3,0} for a 2×2 state *)
DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

(* Return all integer IDs that decode to L×L (perfect-square) arrays
   with ID ≤ maxId.  Identical to vmmc_2d.wl. *)
ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2+1]-1, maxId]]];
      L++];
    ids]

numBeta = 1
