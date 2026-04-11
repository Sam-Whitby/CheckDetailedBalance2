(* ================================================================
   1D weighted-jump dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   State:   flat array of length L.
            State[[i]] = 0 (empty) or k \[Element] {1,...,N} (labeled particle).
            L and N are inferred from the bit string by BitsToState.

   Encoding: bit string \[Rule] integer ID \[Rule] unique labeled array.
             Same bijective encoding as kawasaki_1d.wl / jump_1d_edit.wl.

   Move:    Pick a particle uniformly; draw a displacement d from the
            discrete normal distribution:
              P(d) \[Proportional] w(d) = [Erf((sd+1/2)/(sqrt(2) sigma)) -
                                   Erf((sd-1/2)/(sqrt(2) sigma))] / 2
            where sd is the signed displacement on the ring
            (sd = d if 2d \[LessEqual] L, else sd = d - L).
            This is done via RandomChoice[weights -> displacements],
            which the checker decomposes into sequential Bernoulli trials
            using independent interval-tracking RandomReal[] variables.

            sigma = jumpSigma (set at the bottom of this file).

            The proposal is symmetric: w(d) = w(L-d) by the antisymmetry
            of Erf, so the Metropolis acceptance ratio reduces to
            exp(-beta*dE) exactly, and the algorithm satisfies detailed
            balance.

            If the target site is occupied the move is hard-rejected
            (before the Metropolis step).

   Energy:  Nearest-neighbour pairwise coupling J_ab between particle
            types a < b. Holes (type 0) do not contribute to energy.

   Key difference from jump_1d_edit.wl:
            jump_1d_edit.wl uses rejection sampling (RandomInteger + thinning
            via RandomReal[]) to draw from the normal bin weights.
            This file instead draws directly from the normalised weights via
            RandomChoice[weights -> displacements], which requires the
            checker's Approach-3 weighted-RandomChoice support.  Both
            implementations are mathematically equivalent and both pass
            detailed balance.
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

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Normal bin weight --------------------------------------------------- *)

(* Probability that Normal(0, sigma^2) falls in the unit bin centred on the
   signed displacement sd corresponding to index k on a ring of size L.
   Signed displacement: sd = k if 2k <= L, else sd = k - L.
   w(k) = w(L-k) by the antisymmetry of Erf, so the proposal is symmetric. *)
$normalBinW[k_Integer, sigma_, L_Integer] :=
  With[{sd = If[2 k <= L, k, k - L]},
    (Erf[(sd + 1/2) / (Sqrt[2] sigma)] - Erf[(sd - 1/2) / (Sqrt[2] sigma)]) / 2]


(* ---- Algorithm ----------------------------------------------------------- *)

(* Weighted-choice implementation of a discretised Normal jump.
   RandomChoice[weights -> displacements] draws d directly from the
   normalised bin weights.  The checker (Approach 3) decomposes this into
   n-1 sequential Bernoulli trials, each using a fresh interval-tracking
   RandomReal[] variable, correctly accounting for all conditional
   probabilities.

   The proposal is symmetric (w(d) = w(L-d)), so the forward and reverse
   transition probabilities contributed by the displacement draw cancel
   in the detailed-balance ratio, leaving only the standard Metropolis
   condition exp(-beta*dE). *)
Algorithm[state_List] :=
  Module[{L, particlePositions, nParticles, pickedIdx, pickedPos,
          displacements, weights, d, newPos, newState, dE},
    L = Length[state];
    particlePositions = Flatten[Position[state, _?(# != 0 &)]];
    nParticles = Length[particlePositions];
    If[nParticles == 0, Return[state]];

    (* Choose a particle uniformly at random *)
    pickedIdx = RandomInteger[{1, nParticles}];
    pickedPos = particlePositions[[pickedIdx]];

    (* Draw displacement d in {0,...,L-1} from the discrete normal distribution
       via RandomChoice with normalised bin weights.
       The checker decomposes this into sequential Bernoulli trials. *)
    displacements = Range[0, L - 1];
    weights       = Table[$normalBinW[d, jumpSigma, L], {d, 0, L - 1}];
    d = RandomChoice[weights -> displacements];

    (* Target site under periodic boundary conditions (1-indexed) *)
    newPos = Mod[pickedPos - 1 + d, L] + 1;

    (* Hard rejection if target is occupied (catches d=0 automatically) *)
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

(* Standard deviation of the jump distribution in lattice units.
   Kept as a free symbol during the symbolic check so that the transition
   matrix entries display as exact Erf expressions.  Assigned the value 2
   (in lattice units) during numerical MCMC and animations via the
   "fixedParams" mechanism in symParams. *)
symParams = <|"fixedParams" -> <|jumpSigma -> 2|>|>
