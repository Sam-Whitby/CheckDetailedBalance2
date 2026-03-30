(* ================================================================
   1D uniform-jump dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   State:   flat array of length L.
            State[[i]] = 0 (empty) or k \[Element] {1,...,N} (labeled particle).
            L and N are inferred from the bit string by BitsToState.

   Encoding: same bijective integer encoding as kawasaki_1d.wl /
             jump_1d_edit.wl.

   Move:    Pick a particle uniformly at random.
            Draw a continuous displacement X ~ Uniform[0, L).
            Compute the new site index as:
              newPos = Mod[pickedPos - 1 + X, L] + 1
            which is equivalent to
              Floor[Mod[pickedPos - 1 + X, L]] + 1  (0-indexed floor + 1-offset)
            The checker represents X as a $dbc$contToken["Uniform",0,L,sb].
            The shift (pickedPos-1) propagates through the Plus UpValue,
            and Mod[..., L] folds the range back to [0,L) via the Mod UpValue.
            Floor[...] then collapses the token to a seqBernoulli call with
            uniform weights {1/L, ..., 1/L} over sites {0,...,L-1}, which
            the BFS expands into L branches each with probability 1/L.
            All arithmetic is exact; no floating-point approximation is used.

            If the target site is occupied the move is hard-rejected.

   Energy:  Nearest-neighbour pairwise coupling J_ab between particle
            types a < b. Holes (type 0) do not contribute.

   Why detailed balance holds:
            The displacement is drawn from Uniform[0,L), and every integer
            bin [k, k+1) has probability 1/L, independently of the source
            site.  Thus the proposal probability for a jump of k steps is
            the same in both the forward and reverse directions (since k and
            L-k each have probability 1/L, and the reverse jump has the
            same size modulo L).  The Metropolis factor then ensures that
            the stationary distribution is the Boltzmann distribution.
   ================================================================ *)


(* ---- Bijective integer encoding ----------------------------------------- *)
(* Identical to the encoding in kawasaki_1d.wl / jump_1d_edit.wl.
   Maps every labeled-particle array of any length to a unique integer. *)

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


(* ---- Coupling constants -------------------------------------------------- *)

$pairJ[a_, b_] :=
  If[a == 0 || b == 0, 0,
     ToExpression["J" <> ToString[Min[a,b]] <> ToString[Max[a,b]]]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy -------------------------------------------------------------- *)

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[$pairJ[state[[i]], state[[Mod[i, L] + 1]]], {i, L}]]]


(* ---- Algorithm ----------------------------------------------------------- *)

(* Continuous uniform jump: draw X ~ Uniform[0, L), compute
     newPos = Floor[Mod[pickedPos - 1 + X, L]] + 1
   The checker intercepts RandomReal[{0,L}] as a $dbc$contToken, propagates
   the integer offset (pickedPos-1) through Plus, applies Mod to fold back
   into [0,L), then Floor to collapse to a seqBernoulli with weights {1/L,...}.
   All steps are exact; no floating-point is introduced. *)
Algorithm[state_List] :=
  Module[{L, particlePositions, nParticles, pickedIdx, pickedPos,
          X, newPos, newState, dE},
    L = Length[state];
    particlePositions = Flatten[Position[state, _?(# != 0 &)]];
    nParticles = Length[particlePositions];
    If[nParticles == 0, Return[state]];

    (* Choose a particle uniformly at random *)
    pickedIdx = RandomInteger[{1, nParticles}];
    pickedPos = particlePositions[[pickedIdx]];

    (* Draw a continuous displacement X ~ Uniform[0, L) *)
    X = RandomReal[{0, L}];

    (* New site: shift by X, wrap, discretise (1-indexed) *)
    newPos = Floor[Mod[pickedPos - 1 + X, L]] + 1;

    (* Hard rejection if target is occupied *)
    If[state[[newPos]] != 0, Return[state]];

    (* Build candidate state *)
    newState = ReplacePart[state, {pickedPos -> 0, newPos -> state[[pickedPos]]}];

    (* Metropolis acceptance *)
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
