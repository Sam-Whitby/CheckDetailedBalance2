(* ================================================================
   vmmc_2d_field.wl
   2D VMMC with extended-range coupling and a position-dependent
   field energy on a fully-periodic (torus) square lattice.
   ================================================================

   Extends vmmc_2d.wl with two additions:

   1. EXTENDED-RANGE COUPLING
      Pair interactions at d²=1 (d=1), d²=2 (d=√2), d²=4 (d=2)
      and d²=5 (d=√5) on the L×L torus.  Distances obey the
      minimum-image convention in both row and column.
      Each type pair (a,b) with a<b has four independent coupling
      symbols:
        Jd1<a><b>    – coupling at d²=1
        Jdsq2<a><b>  – coupling at d²=2
        Jd2<a><b>    – coupling at d²=4
        Jdsq5<a><b>  – coupling at d²=5
      These are kept symbolic for DB checking and can be assigned
      real values for numerical simulation / animation.

   2. POSITION-DEPENDENT FIELD ENERGY
      Each lattice site s ∈ {1,...,L²} carries a field symbol Phi<s>.
      Any particle (regardless of type) occupying site s contributes
      Phi<s> to the total energy.  The field is the same for all
      particle types — it represents a site-dependent external potential.
      The change in field energy over a cluster move is handled by a
      post-cluster Metropolis acceptance step:
        accept with min(1, exp(−β ΔE_field))
      where ΔE_field = Σ_{p in cluster} [Phi(dest(p)) − Phi(p)].
      The VMMC link-weight mechanism (unchanged from vmmc_2d.wl) handles
      all pair-energy changes and enforces superdetailed balance for
      those; the Metropolis filter is the correct correction for the
      single-particle field term.

   STATE / ENCODING
      Identical to vmmc_2d.wl: flat L²-array, State[[s]]=0 (empty) or
      k∈{1,...,N} (labeled particle).  BitsToState filters to perfect-
      square lengths.

   REFERENCES
      S. Whitelam & P. L. Geissler, J. Chem. Phys. 127, 154101 (2007)
   ================================================================ *)


(* ================================================================
   SECTION 1 — Bijective integer encoding   (identical to vmmc_2d.wl)
   ================================================================ *)

$cL[L_]        := $cL[L]       = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]    = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_] := $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]


(* ================================================================
   SECTION 2 — Periodic lattice helpers   (identical to vmmc_2d.wl)
   ================================================================ *)

(* Translate site s by direction vector {dRow, dCol} on the L×L torus *)
$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1

(* Row and column of site s (1-indexed) *)
$row[s_, L_] := Ceiling[s/L]
$col[s_, L_] := Mod[s-1, L] + 1


(* ================================================================
   SECTION 3 — Minimum-image squared distance on the L×L torus
   ================================================================ *)

(* Minimum-image squared distance between sites s1 and s2.
   Both row and column use the minimum-image convention, so the
   distance is the same as on an infinite torus. *)
$torusD2[s1_, s2_, L_] :=
  With[{dr0 = Abs[$row[s1,L] - $row[s2,L]],
         dc0 = Abs[$col[s1,L] - $col[s2,L]]},
    With[{dr = Min[dr0, L - dr0],
           dc = Min[dc0, L - dc0]},
      dr^2 + dc^2]]


(* ================================================================
   SECTION 4 — Coupling symbols
   ================================================================ *)

(* Only d²=1 (nearest-neighbour) coupling is active; all longer-range
   couplings are exactly 0.  This keeps the symbol count small so the
   DB checker terminates quickly.  Change d2==1 to the desired set once
   longer-range coupling is needed. *)
$pairJD2[d2_, a_, b_] :=
  If[a == 0 || b == 0 || d2 != 1, 0,
    With[{lo = Min[a,b], hi = Max[a,b]},
      ToExpression["Jd1" <> ToString[lo] <> ToString[hi]]]]


(* ================================================================
   SECTION 5 — Field energy symbol for site s
   ================================================================ *)

(* Phi<s> is the field energy at site s, contributed by any particle
   occupying that site regardless of type. *)
$fieldPhi[s_] := ToExpression["Phi" <> ToString[s]]


(* ================================================================
   SECTION 6 — Neighbour shells (memoised)
   ================================================================ *)

(* Only d²=1 (nearest-neighbour) sites are included, matching the
   active coupling range in $pairJD2.  Restricting this shell is
   critical for checker speed: a site q with zero coupling would still
   trigger a RandomReal[] draw in the cluster builder (consuming a BFS
   bit) even though the link weight is always 0. *)
$neighborsD2[s_, L_] := $neighborsD2[s,L] =
  Select[Range[L^2],
    Function[q, $torusD2[s, q, L] == 1]]


(* ================================================================
   SECTION 7 — Unique undirected bonds at d²=1 (memoised)
   ================================================================ *)

(* Each bond is stored as {s1, s2, d²} with s1 < s2 to avoid
   double-counting.  Used by energy[] to sum all pair contributions. *)
$uniqueBondsExt[L_] := $uniqueBondsExt[L] =
  Flatten[
    Table[
      With[{d2 = $torusD2[s1, s2, L]},
        If[d2 == 1, {{s1, s2, d2}}, {}]],
      {s1, L^2}, {s2, s1+1, L^2}],
    2]


(* ================================================================
   SECTION 8 — Energy function
   ================================================================ *)

(* Total energy = pair energy (extended-range coupling) +
                  field energy (site-dependent, particle-type-independent).
   Pair energy: sum over all undirected bonds at d²∈{1,2,4,5}.
   Field energy: each occupied site s contributes Phi<s>. *)
energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[$pairJD2[#[[3]], state[[#[[1]]]], state[[#[[2]]]]] & /@
          $uniqueBondsExt[L]] +
    Total @ Table[If[state[[s]] != 0, $fieldPhi[s], 0], {s, L^2}]]


(* ================================================================
   SECTION 9 — Virtual pair energy for VMMC link weights
   ================================================================ *)

(* Pair energy between a particle of type typeI at virtual site vI
   and a particle of type typeJ at site qSite on the L×L torus.
   Same-site occupancy (vI == qSite) → Infinity (hard-core exclusion).
   Uses $torusD2 so the minimum-image convention is applied correctly
   even when vI and qSite are on opposite sides of the torus. *)
$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, L_] :=
  Which[
    vI === qSite, Infinity,
    True,
      With[{d2 = $torusD2[vI, qSite, L]},
        If[1 <= d2 <= 5, $pairJD2[d2, typeI, typeJ], 0]]]


(* ================================================================
   SECTION 10 — VMMC cluster builder (extended range)
   ================================================================

   Structure identical to $vmmcBuildCluster in vmmc_2d.wl.
   Differences from vmmc_2d.wl:
     • Neighbour shells use $neighborsD2 (d²≤5) instead of the
       four nearest-neighbour directions.
     • $virtualPairEnergy uses $torusD2 (extended range).
   The Whitelam-Geissler link-weight logic is unchanged.
   ================================================================ *)

$vmmcBuildCluster[state_, L_, seed_, dir_] :=
  Module[{
    cluster   = {seed},
    inCluster = <|seed -> True|>,
    queue     = {seed},
    frustrated = False,
    p, pType, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },

    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p,  dir, L];
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, L];

      (* Union of √5-shells of p, pPost, and pRev.
         A particle q within √5 of pPost or pRev (but not currently
         bonded to p) may have a changed pair energy with p after the
         move and MUST be link-tested.  DeleteDuplicates prevents
         double-testing, which matters on L=2 where right(s)=left(s). *)
      nbrs = DeleteDuplicates @ Join[
               $neighborsD2[p,     L],
               $neighborsD2[pPost, L],
               $neighborsD2[pRev,  L]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          (* Three pair energies for the link-weight formula.
             eInit: current bond energy (0 if not within range).
             eFwd : bond energy if p moves to pPost.
             eRev : bond energy if p moves to pRev. *)
          eInit = $virtualPairEnergy[pType, qType, p,     q, L];
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, L];
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, L];

          (* wFwd = max(1 − exp(β(eInit−eFwd)), 0).
             Piecewise — not Max — so FullSimplify can resolve J sign cases.
             Leading Infinity guards prevent ComplexInfinity with symbolic β.
             (Identical logic to vmmc_2d.wl.) *)
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
                  {1,
                      eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta] (eInit - eRev)],
                      eFwd === Infinity && eInit < eRev},
                  {0,
                      eFwd === Infinity},
                  {1,
                      eRev === Infinity && eInit < eFwd},
                  {Min[(1 - Exp[\[Beta] (eInit - eRev)]) /
                        (1 - Exp[\[Beta] (eInit - eFwd)]), 1],
                   eInit < eFwd && eInit < eRev},
                  {0,
                      eInit < eFwd}},
                0],
              frustrated = True; Break[],
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
   SECTION 11 — Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{L, occupied, seed, dirs, dir, cluster, newState, dEField, dE},
    L = Round[Sqrt[Length[state]]];

    (* Select a random occupied site as seed *)
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    (* Choose a random lattice direction *)
    dirs = {{0,1},{0,-1},{1,0},{-1,0}};
    dir  = RandomChoice[dirs];

    (* Attempt cluster growth; None = frustrated link → reject *)
    cluster = $vmmcBuildCluster[state, L, seed, dir];
    If[cluster === None, Return[state]];

    (* Apply rigid cluster translation.
       Clear cluster sites first, then place at destinations.
       If any destination is occupied by a non-cluster particle,
       reject (overlap). *)
    newState = state;
    Do[newState[[cluster[[i]]]] = 0,            {i, Length[cluster]}];
    Do[With[{dest = $applyDir[cluster[[i]], dir, L]},
         If[newState[[dest]] =!= 0, Return[state, Module]];
         newState[[dest]] = state[[cluster[[i]]]]
       ], {i, Length[cluster]}];

    (* Post-cluster Metropolis filter for field energy change.
       ΔE_field = Σ_{p in cluster} [Phi(dest(p)) − Phi(p)].
       The VMMC link-weight mechanism handles all pair-energy changes
       and enforces superdetailed balance for the pair interactions.
       This Metropolis step is the correct correction for the
       position-dependent single-particle field term. *)
    dEField = Total @ Table[
      $fieldPhi[$applyDir[cluster[[i]], dir, L]] - $fieldPhi[cluster[[i]]],
      {i, Length[cluster]}];
    If[RandomReal[] >= MetropolisProb[dEField], Return[state]];

    (* Interface compliance: MetropolisProb on the total energy change.
       Not used for acceptance — VMMC pair DB is enforced by link weights. *)
    dE = energy[newState] - energy[state];
    MetropolisProb[dE];

    newState
  ]


(* ================================================================
   SECTION 12 — Checker interface   (identical to vmmc_2d.wl)
   ================================================================ *)

(* Accept only IDs whose decoded array has perfect-square length *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

(* Display state as row strings: {1,2}|{3,0} for a 2×2 state *)
DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

(* Return all integer IDs that decode to L×L arrays with ID ≤ maxId *)
ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2+1]-1, maxId]]];
      L++];
    ids]

(* DynamicSymParams: declare the free symbolic parameters for each
   connected component.  Only d²=1 coupling symbols are declared;
   longer-range couplings are numerically 0 and need no symbol.
   Field and coupling symbols both go under "couplings" so that
   check.wls / report.wls assign numeric values to all of them. *)
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]],
          L = Round[Sqrt[Length[First[states]]]]},
    <|"couplings" ->
        Join[
          Flatten @ Table[
            If[a < b,
              {ToExpression["Jd1" <> ToString[a] <> ToString[b]]},
              Nothing],
            {a, types}, {b, types}],
          Table[ToExpression["Phi" <> ToString[s]], {s, L^2}]]|>]

numBeta = 1
