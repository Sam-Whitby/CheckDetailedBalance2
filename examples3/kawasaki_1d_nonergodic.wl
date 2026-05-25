(* ================================================================
   1D Kawasaki dynamics (type-1 only) -- PASSES detailed balance,
   FAILS ergodicity
   ================================================================

   Identical to kawasaki_1d.wl except the swap is only executed if
   exactly one of the two sites holds a type-1 particle and the other
   is a hole (type 0).  Particles with labels ≥ 2 are permanently
   frozen.

   Detailed balance: holds because
     (a) swaps where neither site is type-1 are rejected (T=0 both
         ways), so the DB condition is 0=0 for all such pairs;
     (b) type-1 swaps use standard Metropolis, which satisfies DB.

   Ergodicity: FAILS.  Starting from a state containing any type-2+
   particle (e.g. {1,2,0,0}), those particles never move; the BFS
   reaches only states reachable by type-1 motion, which is a strict
   subset of all (L)_N labeled arrangements.

   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_1d.wl) ----------- *)

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


(* ---- Coupling function (identical to kawasaki_1d.wl) -------------------- *)

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy (identical to kawasaki_1d.wl) -------------------------------- *)

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Algorithm ----------------------------------------------------------- *)

Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L  = Length[state];
    b  = RandomInteger[{0, L - 1}];
    s1 = b + 1; s2 = Mod[b + 1, L] + 1;
    (* Only proceed if exactly one site is type-1 and the other is a hole *)
    If[!((state[[s1]] == 1 && state[[s2]] == 0) ||
         (state[[s2]] == 1 && state[[s1]] == 0)),
      Return[state]];
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
