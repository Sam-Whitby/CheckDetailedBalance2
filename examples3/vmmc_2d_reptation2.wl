(* ================================================================
   2D Reptation VMMC on a periodic square lattice
   ================================================================

   Based on vmmc_2d.wl (Whitelam–Geissler VMMC) with one modification:
   each recruited particle receives an individual move direction rather
   than the same rigid-translation vector as the seed.

   Direction assignment rules for a non-cluster neighbour q of cluster
   particle p (moving from p to pPost = p + pDir):

     (a) Direct collision — q sits exactly at pPost:
           q_dir = pDir   (same as p, conventional VMMC)

     (b) Vacancy-filling — q is a nearest neighbour of p's old site p
         but is NOT pPost:
           q_dir = $dirTo[q, p]   (one step toward p, filling the hole p leaves)

     (c) Shell fallback — q is only in the pPost or pRev neighbour shell:
           q_dir = pDir   (conventional VMMC, needed for correct bond-energy
                           accounting on newly forming / dissolving bonds)

   Rules (a) and (c) reproduce standard VMMC behaviour for those particles.
   Rule (b) produces reptation: neighbours of a vacated site are each
   proposed to fill that hole.  Because all such neighbours share the same
   destination, recruiting more than one of them would cause an intra-cluster
   collision, which aborts the move.

   The VMMC link-probability mechanism (wFwd / wRev / frustrated links) is
   applied identically to vmmc_2d.wl, using each cluster particle's own
   pPost and pRev when evaluating the virtual pair energies.

   Post-cluster Metropolis acceptance:
   In standard VMMC every particle translates rigidly, so intra-cluster
   distances — and hence intra-cluster energies — are unchanged.  Here,
   cluster particles may have different vectors, so intra-cluster bond
   energies can change.  A Metropolis step on dE_intra corrects for this.
   Cluster–noncluster energy changes are still handled entirely by the
   link probabilities (as in standard VMMC).

   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_2d.wl) ----------- *)

$cL[L_]        := $cL[L]     = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]  = Sum[$cL[l], {l, 0, L - 1}]
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


(* ---- 2D lattice helpers -------------------------------------------------- *)

$right2D[s_, L_] := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c,L] + 1]
$down2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r,L]*L + c]
$left2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c-2,L] + 1]
$up2D[s_, L_]    := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r-2,L]*L + c]

$allNeighbors2D[s_, L_] := {$right2D[s,L], $left2D[s,L], $down2D[s,L], $up2D[s,L]}

$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1

(* Unit direction vector that steps from [from] to adjacent site [to].
   Returns {0,0} if [to] is not a nearest neighbour (should never occur
   in normal use, where the caller has already confirmed adjacency). *)
$dirTo[from_, to_, L_] :=
  Which[
    $right2D[from, L] === to, {0,  1},
    $left2D[from, L]  === to, {0, -1},
    $down2D[from, L]  === to, {1,  0},
    $up2D[from, L]    === to, {-1, 0},
    True,                     {0,  0}
  ]


(* ---- Coupling function --------------------------------------------------- *)

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy -------------------------------------------------------------- *)

$uniqueBonds2D[L_] := $uniqueBonds2D[L] =
  DeleteDuplicates[Sort /@ Flatten[
    Table[{{s, $right2D[s,L]}, {s, $down2D[s,L]}}, {s, L^2}], 1]]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[couplingJ[state[[#[[1]]]], state[[#[[2]]]], 1] & /@ $uniqueBonds2D[L]]]


(* ---- VMMC helpers -------------------------------------------------------- *)

$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, L_] :=
  Which[
    vI === qSite,                           Infinity,
    MemberQ[$allNeighbors2D[vI, L], qSite], couplingJ[typeI, typeJ, 1],
    True,                                   0]


(* ---- Reptation VMMC cluster builder ------------------------------------- *)

(* Returns {cluster, clusterDir} on success, or None if any frustrated link
   or intra-cluster collision is detected.

   clusterDir is an Association from each cluster site to its move vector.
   claimedPos tracks the new (destination) sites already assigned to cluster
   members; a second particle trying to move to the same site aborts the move. *)

$vmmcBuildCluster[state_, L_, seed_, dir_] :=
  Module[{
    cluster    = {seed},
    inCluster  = <|seed -> True|>,
    clusterDir = <|seed -> dir|>,
    claimedPos = <|$applyDir[seed, dir, L] -> True|>,
    queue      = {seed},
    frustrated = False,
    p, pType, pDir, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2, qDir, qDest
  },
    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pDir  = clusterDir[p];
      pPost = $applyDir[p, pDir, L];
      pRev  = $applyDir[p, {-pDir[[1]], -pDir[[2]]}, L];
      (* Same combined neighbour shell as vmmc_2d.wl *)
      nbrs  = DeleteDuplicates[Join[
                $allNeighbors2D[p, L],
                $allNeighbors2D[pPost, L],
                $allNeighbors2D[pRev, L]]];
      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];
          (* Assign direction to q based on its geometric relationship to p *)
          qDir = Which[
            q === pPost,
              pDir,                  (* (a) direct collision: same direction *)
            MemberQ[$allNeighbors2D[p, L], q],
              $dirTo[q, p, L],       (* (b) vacancy-filling: step into p's old site *)
            True,
              pDir                   (* (c) shell fallback: conventional VMMC *)
          ];
          (* Link probabilities use p's own pPost / pRev — identical to vmmc_2d.wl *)
          eInit = $virtualPairEnergy[pType, qType, p,     q, L];
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, L];
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, L];
          wFwd = Piecewise[{
              {1,                               eFwd === Infinity},
              {1 - Exp[\[Beta] (eInit - eFwd)], eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                               eRev === Infinity},
              {1 - Exp[\[Beta] (eInit - eRev)], eInit < eRev}},
            0];
          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1,                               eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta] (eInit - eRev)], eFwd === Infinity && eInit < eRev},
                  {0,                               eFwd === Infinity},
                  {1,                               eRev === Infinity && eInit < eFwd},
                  {Min[(1 - Exp[\[Beta] (eInit - eRev)]) / (1 - Exp[\[Beta] (eInit - eFwd)]), 1],
                   eInit < eFwd && eInit < eRev},
                  {0,                               eInit < eFwd}},
                0],
              frustrated = True; Break[],    (* frustrated link: reject *)
              (* Intra-cluster collision check *)
              qDest = $applyDir[q, qDir, L];
              If[KeyExistsQ[claimedPos, qDest],
                frustrated = True; Break[], (* two particles target same site *)
                AppendTo[cluster, q];
                inCluster[q]  = True;
                clusterDir[q] = qDir;
                claimedPos[qDest] = True;
                AppendTo[queue, q]
              ]
            ]
          ]
        ],
        {k, Length[nbrs]}
      ]
    ];
    If[frustrated, None, {cluster, clusterDir}]
  ]


(* ---- Algorithm ----------------------------------------------------------- *)

Algorithm[state_List] :=
  Module[{
    L, occupied, seed, dirs, dir,
    result, cluster, clusterDir, n,
    newState, dest, dEIntra
  },
    L        = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];
    dirs = {{0,1},{0,-1},{1,0},{-1,0}};
    dir  = RandomChoice[dirs];

    result = $vmmcBuildCluster[state, L, seed, dir];
    If[result === None, Return[state]];
    {cluster, clusterDir} = result;

    (* Apply moves: clear all old positions first, then place at destinations *)
    newState = state;
    Do[newState[[cluster[[i]]]] = 0,                              {i, Length[cluster]}];
    Do[dest = $applyDir[cluster[[i]], clusterDir[cluster[[i]]], L];
       If[newState[[dest]] =!= 0, Return[state, Module]];
       newState[[dest]] = state[[cluster[[i]]]],                  {i, Length[cluster]}];

    (* Intra-cluster Metropolis correction.
       The VMMC link probabilities account for all cluster–noncluster energy
       changes (as in standard VMMC).  Intra-cluster energy changes are
       non-zero here because cluster particles do not share a single rigid
       translation, so relative positions within the cluster can change.
       MetropolisProb[dE_intra] corrects for this. *)
    n = Length[cluster];
    dEIntra = Sum[
      With[{si = cluster[[i]], sj = cluster[[j]],
            ti = $applyDir[cluster[[i]], clusterDir[cluster[[i]]], L],
            tj = $applyDir[cluster[[j]], clusterDir[cluster[[j]]], L]},
        If[MemberQ[$allNeighbors2D[ti, L], tj], couplingJ[state[[si]], state[[sj]], 1], 0] -
        If[MemberQ[$allNeighbors2D[si, L], sj], couplingJ[state[[si]], state[[sj]], 1], 0]
      ],
      {i, 1, n}, {j, i+1, n}
    ];

    If[RandomReal[] < MetropolisProb[dEIntra], newState, state]
  ]


(* ---- Checker interface (identical to vmmc_2d.wl) ------------------------ *)

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
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]

numBeta = 1
