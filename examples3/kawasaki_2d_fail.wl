(* ================================================================
   2D Kawasaki with wrong-sign acceptance -- FAILS detailed balance
   ================================================================

   Identical to kawasaki_2d.wl except for one line in Algorithm:

     CORRECT:  dE = energy[newState] - energy[state]
     BUGGY:    dE = energy[state] - energy[newState]

   Same bug as kawasaki_1d_fail.wl: the sign of dE is reversed, so
   the algorithm prefers high-energy states (uphill moves always
   accepted, downhill moves penalised). Detailed balance is violated
   for any pair of states with different energies.
   ================================================================ *)


(* ---- Encoding and lattice helpers (identical to kawasaki_2d.wl) --------- *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

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

$right2D[s_, L_] := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c,L] + 1]
$down2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r,L]*L + c]

$bond2D[L_, b_] :=
  If[b < L*L,
    {b+1, $right2D[b+1, L]},
    With[{s = b - L*L + 1}, {s, $down2D[s, L]}]]

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[Table[
      couplingJ[state[[s]], state[[$right2D[s, L]]], 1] +
      couplingJ[state[[s]], state[[$down2D[s, L]]], 1],
      {s, Length[state]}]]]


(* ---- Buggy algorithm: dE computed with wrong sign ----------------------- *)

Algorithm[state_List] :=
  Module[{L, b, sites, s1, s2, newState, dE},
    L     = Round[Sqrt[Length[state]]];
    b     = RandomInteger[{0, 2*L*L - 1}];
    sites = $bond2D[L, b]; s1 = sites[[1]]; s2 = sites[[2]];
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[state] - energy[newState];  (* BUG: sign reversed -- prefers high energy *)
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

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

numBeta = 1
