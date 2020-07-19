from collections  import namedtuple
from itertools    import count, izip
from glob         import glob
from time         import time, sleep

from synthesis    import Synthesis
from utils        import Fold
from kmerSetDB    import kmerSetDB
from kmerSetArray import kmerSetArray

import sys
import random
import uuid
import math
import textwrap

import utils
import makerchecks
import projector
import berno
import finder

from Bio.SeqUtils import MeltingTemp


class NRPMaker(object):

    def __init__(self, part_type='RNA', seed=None):
        self.iupac_table = {
            'A': {'A'},
            'C': {'C'},
            'G': {'G'},
            'T': {'T'},
            'R': {'A', 'G'},
            'Y': {'C', 'T'},
            'S': {'G', 'C'},
            'W': {'A', 'T'},
            'K': {'G', 'T'},
            'M': {'A', 'C'},
            'B': {'C', 'G', 'T'},
            'V': {'A', 'C', 'G'},
            'D': {'A', 'G', 'T'},
            'H': {'A', 'C', 'T'},
            'N': {'A', 'T', 'G', 'C'}
        }
        self.compl_iupac_map = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'R': 'Y',
            'Y': 'R',
            'S': 'W',
            'W': 'S',
            'K': 'M',
            'M': 'K',
            'B': 'V',
            'V': 'B',
            'D': 'H',
            'H': 'D',
            'N': 'N'
        }
        self.compl_dict = {
            'A': set(['T']),
            'G': set(['C', 'T']),
            'C': set(['G']),
            'T': set(['A', 'G'])
        }
        self.synthesis  = Synthesis()
        self.fold       = Fold(part_type=part_type)        
        self.proj_id    = str(uuid.uuid4())
        self.kmer_db    = None
        self.background = None
        
        if not seed is None and isinstance(seed, int):
            self.seed = seed
        else:
            self.seed = random.random()
        random.seed(self.seed)

    def _get_adjusted_struct(self, struct, seq):
        if struct is None:
            return '.'*len(seq)
        elif len(seq) < len(struct):
            return struct[:len(seq)]
        elif len(seq) > len(struct):
            return ''.join([struct, '.'*(len(seq) - len(struct))])
        return struct

    def _get_meta_struct(self, struct):
        struct_spec   = namedtuple(
            'struct_spec', 'struct paired_dict rev_paired_dict unpaired_set folding inversefolding')
        pairing_stack = []
        meta_struct   = struct_spec(
            struct=struct.replace(
                'x', '.'),
            paired_dict={},
            rev_paired_dict={},
            unpaired_set=set(),
            folding={},
            inversefolding={})

        for index, nt in enumerate(struct):

            if nt == '(':
                pairing_stack.append(index)

            elif nt == ')':

                try:
                    closure = pairing_stack.pop()
                except:
                    raise ValueError(
                        ' [X] Unbalanced parentheses in structure constraint')

                meta_struct.paired_dict[closure]   = index
                meta_struct.rev_paired_dict[index] = closure

            elif nt == 'x':
                meta_struct.unpaired_set.add(index)

        if pairing_stack:
            raise ValueError(
                ' [X] Unbalanced parentheses in structure constraint')

        if meta_struct.paired_dict:
            meta_struct.folding['status'] = True
            meta_struct.inversefolding['status'] = True
        elif meta_struct.unpaired_set:
            meta_struct.folding['status'] = True
            meta_struct.inversefolding['status'] = False
        else:
            meta_struct.folding['status'] = False
            meta_struct.inversefolding['status'] = False

        return meta_struct

    def _get_meta_seq(self, seq, meta_struct):
        seq = list(seq)
        meta_seq = [None] * len(seq)

        # Normalize Meta Sequence
        for i in xrange(len(seq)):            
            try:
                if i in meta_struct.paired_dict:
                    j = meta_struct.paired_dict[i]
                    if len(self.iupac_table[seq[j]]) < len(self.iupac_table[seq[i]]):
                        seq[i] = self.compl_iupac_map[seq[j]]
                meta_seq[i] = set(self.iupac_table[seq[i]])
            except:
                raise ValueError(
                    ' [X] Invalid IUPAC code at index {} in sequence constraint'.format(i))

        return tuple(meta_seq)

    def _reset_candidate_kmer_set(
        self,
        candidate,
        kmer_set,
        i):
        candidate[i] = '-'
        kmer_set[i]  = ' '

    def _freeze_point(
        self,
        candidate,
        tried_set,
        kmer_set,
        i,
        k):
        j = i
        while j > k:
            self._reset_candidate_kmer_set(
                candidate,
                kmer_set,
                j)
            tried_set[j] = set()
            # tried_set[j] = set(candidate[j])
            j -= 1
        else:
            self._reset_candidate_kmer_set(
                candidate,
                kmer_set,
                j)
        i = j
        return i

    # FUNCTION TO TRACEBACK TO OPTIMAL LOCATION UPON CONFLICT
    # Note: Compute these values for all position at the beginning 
    # for the given sequence and structure constraints.
    def _get_rollback_index(
        self,
        meta_seq,
        meta_struct,
        i,
        homology):
        roll_back_index = None
        # See if any position in last homology places has potential for change
        j = i-1
        while j > i - homology and j > -1:
            # Case ( or . without fixed base in sequence constraint
            if not (j in meta_struct.rev_paired_dict) and (len(meta_seq[j]) != 1):
                roll_back_index = j
                break
            j -= 1
        # Else go to the ( for the first ) in the last homology places
        if not roll_back_index:
            j = max(0, i - homology + 1)
            k = None
            while j <= i:
                if j in meta_struct.rev_paired_dict:
                    k = j
                    break
                j += 1
            roll_back_index = meta_struct.rev_paired_dict[k]
        return roll_back_index

    def _roll_back(
        self,
        meta_seq,
        meta_struct,
        i,
        homology,
        candidate,
        tried_set,
        kmer_set,
        rbi=None):
        # Default Case: Explicit roll back
        if not rbi is None:
            roll_back_index = rbi
            self._freeze_point(
                candidate, tried_set, kmer_set, i, k=roll_back_index)
            return rbi

        # Case 1: i is not paired upstream
        if not i in meta_struct.rev_paired_dict:
            self._reset_candidate_kmer_set(
                candidate, kmer_set, i)

        # Case 2: i is paired upstream
        else:
            roll_back_index = self._get_rollback_index(
                meta_seq, meta_struct, i, homology)
            i = self._freeze_point(
                candidate, tried_set, kmer_set, i, k=roll_back_index)
        
        return i

    def _get_local_roll_back_index(
        self,
        candidate,
        i,
        local_model_fn,
        verbose):
        # Prep candidate
        candidate_str = ''.join(candidate[:i+1])
        
        # Try to evaluate the local_model_fn on candidate_str
        outcome = True
        try:
            outcome = local_model_fn(candidate_str)
            if outcome in [True, False]:
                state, index = outcome, i
            else:
                state, index = outcome
        except Exception as e:
            print '\n Local Model fn. failed to evaluate partial path: {}\n'.format(
                candidate_str)
            raise e # No intelligence, halt everything!

        # State satisfactory?
        try:
            assert state in [True, False]
        except Exception as e:
            print '\n Local Model fn. returned a non-boolean evaluation: {}\n'.format(
                state)
            raise e

        # Index satisfactory?
        try:
            if state == False:
                index = int(index)
                assert index <= i
        except Exception as e:
            print '\n Local Model fn. returned a non-integer traceback index: {}\n'.format(
                index)
            raise e

        # No conflict found!
        return None

    def _get_non_coding_candidate(
        self,
        meta_seq,
        meta_struct,
        homology,
        local_model_fn,
        verbose,
        jump=False,
        start_seq=None,
        allow_internal_repeat=False):

        # Setup the data structures
        candidate = ['-'] * len(meta_seq) if not start_seq else list(start_seq)
        tried_set = [set() for _ in xrange(len(meta_seq))]
        kmer_set  = kmerSetArray(size=len(meta_seq))

        # Setup indexing
        i               = 0
        roll_back_count = 0

        # Main backtracking code
        while -1 < i < len(meta_seq):

            # Jumping out of iteration
            if jump and roll_back_count == homology:
                candidate = None
                break

            # Try to build a candidate
            if candidate[i] == '-':

                # Forward phase - A nucleotide may be chosen
                if len(tried_set[i]) < len(meta_seq[i]):
                    # Reset roll_back_count
                    roll_back_count = 0

                    # Case ( and .
                    if not i in meta_struct.rev_paired_dict:
                        candidate[i] = random.choice(
                            list(meta_seq[i]-tried_set[i]))
                        tried_set[i].add(candidate[i])
                        # Case (
                        if i in meta_struct.paired_dict:
                            candidate[meta_struct.paired_dict[i]] = utils.get_revcomp(
                                candidate[i])
                            tried_set[meta_struct.paired_dict[i]] = set(
                                candidate[meta_struct.paired_dict[i]])
                    # Case )
                    else:
                        candidate[i] = utils.get_comp(
                            candidate[meta_struct.rev_paired_dict[i]])
                        tried_set[i] = set(candidate[i])

                # Backward phase - Nucleotide choices exhausted, so traceback
                else:
                    # Update roll_back_count
                    roll_back_count += 1

                    # Reset and roll back
                    tried_set[i] = set()
                    kmer_set[i]  = ' '
                    if i > 0:
                        # Clear stuff at current index
                        self._reset_candidate_kmer_set(
                            candidate,
                            kmer_set,
                            i)
                        # Traceback to previous index
                        # since current index done
                        i = self._roll_back(
                            meta_seq,
                            meta_struct,
                            i-1,
                            homology,
                            candidate,
                            tried_set,
                            kmer_set) + 1
                    i -= 1
                    continue

            # See if built candidate[i] is valid
            # Case ( and .
            elif not i in meta_struct.rev_paired_dict:
                # Wrong nucleotide selected
                if candidate[i] not in meta_seq[i]:
                    candidate[i] = '-'
                    kmer_set[i]  = ' '
                    tried_set[i] = set()
                    continue
            # Case )
            elif i in meta_struct.rev_paired_dict:
                j = meta_struct.rev_paired_dict[i]
                # Paired bases are not complementary
                if utils.get_comp(candidate[i]) != candidate[j]:
                    candidate[i] = '-'
                    kmer_set[i]  = ' '
                    tried_set[i] = set()
                    continue

            # See if the local model function is violated
            if local_model_fn:
                rbi = self._get_local_roll_back_index(
                    candidate, i, local_model_fn, verbose)
                # Model function violated, and
                # a traceback location was determined
                if not rbi is None: 
                    roll_back_count = 0 if rbi < i else roll_back_count
                    i = self._roll_back(
                        meta_seq,
                        meta_struct,
                        i,
                        homology,
                        candidate,
                        tried_set,
                        kmer_set,
                        rbi)
                    continue

            # Handle internal and shared repeats
            if i >= homology-1:
                # Get the kmer/rmer
                kmer = ''.join(candidate[i-homology+1:i+1])
                rmer = utils.get_revcomp(kmer)
                mmer = min(kmer, rmer)

                # Are either of these mers seen previously?
                mmer_seen = False

                # Case: mmer is a shared repeat with
                #       a previous part
                if mmer in self.kmer_db:
                    mmer_seen = True

                # Case: mmer is a shared repeat with
                #       background
                if not self.background is None:
                    if mmer in self.background:
                        mmer_seen = True

                # Case: kmer/rkmer is an internal
                #       repeat to current part
                elif not allow_internal_repeat:
                    # Direct repeat
                    if kmer in kmer_set:
                        mmer_seen =True
                    # Inverted repeat
                    if rmer in kmer_set:
                        mmer_seen=True

                # Traceback to eliminate repeat
                if mmer_seen:
                    i = self._roll_back(
                        meta_seq,
                        meta_struct,
                        i,
                        homology,
                        candidate, 
                        tried_set,
                        kmer_set)
                    continue
                
                else:
                    kmer_set[i] = kmer

            # Roll forward
            i += 1

        # Prepare to return candidate
        if candidate is None or '-' in candidate:
            return None
        else:
            return ''.join(candidate)

    def _get_opt_pass_count(self, meta_struct, synth_opt, global_model_fn):
        opt_criteria_count = 1 # Since candidate must at least be non-repetitive
        if meta_struct.folding['status']:
            opt_criteria_count += 1
        if synth_opt:
            opt_criteria_count += 1
        if global_model_fn:
            opt_criteria_count += 1
        return opt_criteria_count

    # Diagnostic function -- Shouldn't trigger on experimental changes
    def _is_non_coding_construction_verified(self, meta_seq, meta_struct, candidate):
        i = 0
        while i < len(candidate):
            if not candidate[i] in meta_seq[i]:
                return False
            if i in meta_struct.paired_dict:
                j = meta_struct.paired_dict[i]
                if not candidate[j] in self.compl_dict[candidate[i]]:
                    # print '{}th Nucleotide {} not compementary to {}th Nucleotide {}'.format(
                    #     i,
                    #     candidate[i],
                    #     j,
                    #     candidate[j])
                    return False
            i += 1
        return True

    def _is_synthesis_verified(self, candidate):
        return self.synthesis.evaluate(candidate)

    def _is_structure_verified(self, meta_struct, struct_type, candidate):
        struct_satisfied = 0

        # Decide which structure criteria to fulfill
        if struct_type == 'centroid':
            candidate_structs = [self.fold.evaluate_centroid(
                candidate)]
        elif struct_type == 'both':
            candidate_structs = [self.fold.evaluate_mfe(
                candidate)]
            candidate_structs.append(self.fold.evaluate_centroid(
                candidate))
        else:
            candidate_structs = [self.fold.evaluate_mfe(
                candidate)]

        for candidate_struct in candidate_structs:
            vienna_meta_struct = self._get_meta_struct(
                candidate_struct)
            # Ensure all forbidden base indices unpaired
            for bp_closed in meta_struct.unpaired_set:
                if not vienna_meta_struct.struct[bp_closed] == '.':
                    return False
            # Ensure all paired base indices paired as desired
            for bp_open in meta_struct.paired_dict:
                if not bp_open in vienna_meta_struct.paired_dict:
                    break
                else:
                    if meta_struct.paired_dict[bp_open] != vienna_meta_struct.paired_dict[bp_open]:
                        break
            else:
                struct_satisfied += 1

        # All structure criteria satisfied
        if struct_satisfied == len(candidate_structs):
            return True
        else:
            return False

    def _get_inverse_fold_candidate(self, start_seq, meta_seq, meta_struct):
        inverse_fold_seq = ''.join(char.lower() if len(meta_seq[i]) == 1 else char for i,char in enumerate(start_seq))
        return self.fold.design(
            seq=inverse_fold_seq, struct=meta_struct.struct).upper()

    def _get_verified_non_coding_candidate(self, 
        homology,
        meta_seq,
        meta_struct,
        struct_type,
        synth_opt,
        local_model_fn,
        global_model_fn,
        jump_count,
        fail_count,
        verbose,
        abortion,
        allow_internal_repeat=False):

        # Setup counts and variables
        current_jump_count = 0
        current_fail_count = 0
        struct_fail_count  = 0
        synth_fail_count   = 0
        model_fail_count   = 0
        seed_seq           = None
        verified_candidate = None
        opt_pass_count     = self._get_opt_pass_count(
            meta_struct,
            synth_opt,
            global_model_fn)


        while not verified_candidate:

            # Try to get a non-repetitive candidate
            start_seq = self._get_non_coding_candidate(
                meta_seq=meta_seq,
                meta_struct=meta_struct,
                homology=homology,
                local_model_fn=local_model_fn,
                verbose=verbose,
                jump=current_jump_count < jump_count,
                start_seq=seed_seq,
                allow_internal_repeat=allow_internal_repeat)
            candidate = start_seq
            if start_seq:
                # If structure unmatched but constraint has base pairings, go for inverse-repair strategy
                if meta_struct.inversefolding['status'] and not self._is_structure_verified(
                    meta_struct, struct_type, start_seq):
                    inverse_fold_seq = self._get_inverse_fold_candidate(
                        start_seq,
                        meta_seq,
                        meta_struct)
                    candidate = self._get_non_coding_candidate(
                        meta_seq=meta_seq,
                        meta_struct=meta_struct,
                        homology=homology,
                        local_model_fn=local_model_fn,
                        verbose=verbose,
                        jump=current_jump_count < jump_count,
                        start_seq=inverse_fold_seq,
                        allow_internal_repeat=allow_internal_repeat)

            # print candidate
            opt_count = 0

            # If valid candidate then process
            if candidate:
                current_jump_count = 0
                opt_count += 1

                # Diagnostic block -- Shouldn't trigger on experimental changes
                if not self._is_non_coding_construction_verified(
                    meta_seq,
                    meta_struct,
                    candidate):
                    current_fail_count += 1
                    raise Exception(
                        'Part construction incorrect due to code/logic changes! Report issue.')
                else:
                    pass

                # Synthesis optimization            
                if synth_opt:
                    if self._is_synthesis_verified(candidate):
                        opt_count += 1
                    else:
                        synth_fail_count += 1
                
                # Global model optimization
                if global_model_fn:

                    # Try to evaluate the global_model_fn on candidate_str
                    outcome = True
                    try:
                        outcome = global_model_fn(candidate)
                    except Exception as e:
                        print '\n Global Model fn. failed to evaluate complete path: {}\n'.format(
                            candidate)
                        raise e

                    # Outcome satisfactory?
                    try:
                        assert outcome in [True, False]
                    except Execution as e:
                        print '\n Global Model fn. returned a non-boolean state: {}\n'.format(
                            outcome)
                        raise e

                    # Process outcome
                    if outcome: # True
                        opt_count += 1
                    else:       # False
                        model_fail_count += 1

                # Structural optimization
                if meta_struct.folding['status']:
                    if self._is_structure_verified(
                        meta_struct,
                        struct_type,
                        candidate):
                        opt_count += 1
                    else:
                        struct_fail_count += 1

                # Did everything get optimized?
                if opt_count == opt_pass_count:
                    verified_candidate = candidate
                else:
                    current_fail_count += 1

                # Failure count exceeded, terminate
                if current_fail_count == fail_count:
                    break

            # No candidate produced
            else:
                # No jumps made yet no non-repetitive candidate found
                if current_jump_count >= jump_count:
                    break
                # Increase current_jump_count 
                else:
                    current_jump_count += 1
                # Abortion limit reached?
                if abortion and current_jump_count >= jump_count:
                    break

        if int(verbose) > 1:
            print '\n  [seq fails] {}, [struct fails] {}, [synth fails] {}, [global fails] {}, [opt fails] {}'.format(
                current_jump_count,
                struct_fail_count,
                synth_fail_count,
                model_fail_count,
                current_fail_count)

        return verified_candidate, current_jump_count+1, current_fail_count+1

    def _get_non_coding_nrps(
        self,
        homology,
        seq,
        struct,
        struct_type,
        target,
        synth_opt,
        local_model_fn,
        global_model_fn,
        jump_count,
        fail_count,
        verbose,
        abortion,
        allow_internal_repeat=False):

        # Setup structures
        seq = seq.upper().replace('U', 'T')
        meta_struct = self._get_meta_struct(
            struct=self._get_adjusted_struct(struct, seq))
        meta_seq = self._get_meta_seq(
            seq=seq.upper().replace('U', 'T'),
            meta_struct=meta_struct)

        seq_count  = 0
        iter_count = 0
        time_sum   = 0.0
        begin_time = time()
        break_flag = False

        # Setup Bernoulli Success model
        total_jump_trials    = jump_count
        total_jump_successes = 1
        curr_jump_prob       = berno.get_prob(
            trials=total_jump_trials,
            success=total_jump_successes)
        curr_jump_trial      = jump_count
        total_fail_trials    = fail_count
        total_fail_successes = 1
        curr_fail_prob       = berno.get_prob(
            trials=total_fail_trials,
            success=total_fail_successes)
        curr_fail_trial      = fail_count

        # Stream parts until completion
        while True:
            
            t0 = time()

            candidate, curr_jump_trial, curr_fail_trial = self._get_verified_non_coding_candidate(
                homology,
                meta_seq,
                meta_struct,
                struct_type,
                synth_opt,
                local_model_fn,
                global_model_fn,
                curr_jump_trial,
                curr_fail_trial,
                verbose,
                abortion,
                allow_internal_repeat)

            if candidate is None:
                break_flag = True

            # Got a canidate -- will try again
            if not break_flag:
                final_candidate = candidate
                update_status   = True
                try:
                    self.kmer_db.add(candidate)
                except Exception as E:
                    update_status = False

                if not update_status: # Memory full
                    if verbose:
                        print '[ERROR] Memory Full ... Breaking Loop'
                        yield update_status

                seq_count  += 1
                time_sum   += time()-t0
                iter_count += 1

                if verbose:
                    print ' [part] {}, [{}-mers] {}, [iter {}] {:.2f}s, [avg] {:.2f}s, [total time] {:.2f}h'.format(
                        seq_count,
                        homology,
                        len(self.kmer_db),
                        iter_count,
                        time()-t0, time_sum / iter_count,
                        (time() - begin_time) / 3600.0)
                
                yield final_candidate

                # No more parts required
                if seq_count == target:
                    break
            # No more candidates to build
            else:
                yield candidate
                break

            # Update failure limits based on Bernoulli Success model
            total_jump_trials    += curr_jump_trial
            total_jump_successes += 1
            curr_jump_prob       =  berno.get_prob(trials=total_jump_trials, success=total_jump_successes)
            curr_jump_trial      =  berno.get_trials(prob=curr_jump_prob)
            total_fail_trials    += curr_fail_trial
            total_fail_successes += 1
            curr_fail_prob       =  berno.get_prob(trials=total_fail_trials, success=total_fail_successes)
            curr_fail_trial      =  berno.get_trials(prob=curr_fail_prob)

    def _check_maker_constraints(
        self,
        seq,
        struct,
        allow_internal_repeat,
        target,
        homology):
        # Sequence Legality 1
        if not isinstance(seq, str):
            print '\n [ERROR]    Sequence Constraint must be a string, not \'{}\''.format(seq)
            print ' [SOLUTION] Try correcting Sequence Constraint\n'
            return False
        # Sequence Legality 2
        seq_legal, seq_illegal_chars = makerchecks.is_seq_constr_legal(seq)
        if not seq_legal:
            print '\n [ERROR]    Sequence Constraint is not legal due to chars: {}'.format(seq_illegal_chars)
            print ' [SOLUTION] Try correcting Sequence Constraint\n'
            return False
        # Structure Legality 1
        if not isinstance(struct, str):
            print '\n [ERROR]    Structure Constraint must be a string, not \'{}\''.format(struct)
            print ' [SOLUTION] Try correcting Structure Constraint\n'
            return False
        # Lmax Legality 1
        if not isinstance(homology, int):
            print '\n [ERROR]    Lmax must be an integer, not \'{}\''.format(homology-1)
            print ' [SOLUTION] Try correcting Lmax\n'
            return False
        # Lmax Legality 2
        if homology-1 < 5:
            print '\n [ERROR]    Lmax must be greater than 4, not \'{}\''.format(homology-1)
            print ' [SOLUTION] Try correcting Lmax\n'
            return False
        # Target Size Legality 1
        if not isinstance(target, int):
            print '\n [ERROR]    Target Size must be an integer, not \'{}\''.format(target)
            print ' [SOLUTION] Try correcting Target Size\n'
            return False
        # Target Size Legality 2
        if target < 1:
            print '\n [ERROR]    Target Size must be greater than 0, not \'{}\''.format(target)
            print ' [SOLUTION] Try correcting Target Size\n'
            return False
        # Sequence Sufficiency 1 -- soft check
        seq_sufficient, constrained_motif_locs = makerchecks.is_seq_constr_sufficient(seq, homology, target)
        if not seq_sufficient:
            print '\n [WARNING]  Target size of {} is unreachable from given Sequence Constraint and Lmax of {}'.format(target, homology-1)
            print ' [WARNING]  >> Overly constrained motifs at locations: {}'.format(constrained_motif_locs)
            print ' [WARNING]  Fewer parts will be generated'
        # Structure Legality 2
        struct_legal, unclosed, unopened, invalid = makerchecks.is_structure_valid(struct)
        if not struct_legal:
            print '\n [ERROR]    Structure Constraint is illegal or unbalanced'
            if unclosed:
                print ' [ERROR]    >> Unclosed bases at locations: {}'.format(unclosed)
            if unopened:
                print ' [ERROR]    >> Unopened bases at locations: {}'.format(unopened)
            if invalid:
                print ' [ERROR]    >> Invalid characters at locations: {}'.format(invalid)
            print ' [SOLUTION] Try correcting Structure Constraint\n'
            return False
        # Internal Repeats Legality
        if not allow_internal_repeat in [True, False]:
            print '\n [ERROR]    Internal Repeat must be boolean, not \'{}\''.format(allow_internal_repeat)
            print ' [SOLUTION] Try correcting Internal Repeat\n'
            return False
        # Structure Sufficiency 2
        if not allow_internal_repeat:
            struct_sufficient, long_hairpins = makerchecks.is_structure_not_conflict(struct, homology)
            if not struct_sufficient:
                print '\n [ERROR] Structure Constraint is insufficient based on given Lmax'
                for long_hairpin in long_hairpins:
                    print ' [ERROR] >> Long hairpin between: {}'.format(long_hairpin)
                    print ' [SOLUTION] Try relaxing Structure Constraint or setting allow_internal_repeat=True\n'
                return False
        return True

    def _check_maker_inputs(
        self,
        part_type,
        struct_type,
        synth_opt,
        jump_count,
        fail_count,
        output_file,
        verbose):
        # Part Type Legality
        if not part_type in ['DNA', 'RNA']:
            print '\n [ERROR]    Part Type must be \'RNA\' or \'DNA\', not \'{}\''.format(part_type)
            print ' [SOLUTION] Try correcting Part Type\n'
            return False
        # Struct Type Legality
        if not struct_type in ['mfe', 'centroid', 'both']:
            print '\n [ERROR]    Struct Type must be \'mfe\', \'centroid\' or \'both\', not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Part Type\n'
            return False
        # Synth Optimization Legality
        if not synth_opt in [True, False]:
            print '\n [ERROR]    Synth Opt must be True or False, not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Synth Opt\n'
            return False
        # Jump Count Legality 1
        if not isinstance(jump_count, int):
            print '\n [ERROR]    Jump Count must be an integer, not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Jump Count\n'
            return False
        # Jump Count Legality 2
        if jump_count < 0:
            print '\n [ERROR]    Jump Count must be an positive, not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Jump Count\n'
            return False
        # Fail Count Legality 1
        if not isinstance(fail_count, int):
            print '\n [ERROR]    Fail Count must be an integer, not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Fail Count\n'
            return False
        # Fail Count Legality 2
        if fail_count < 0:
            print '\n [ERROR]    Fail Count must be an positive, not \'{}\''.format(struct_type)
            print ' [SOLUTION] Try correcting Fail Count\n'
            return False
        # Output File Legality
        if not output_file is None:
            if not isinstance(output_file, str):
                print '\n [ERROR]    Output File must be a string or None, not \'{}\''.format(struct_type)
                print ' [SOLUTION] Try correcting Output File\n'
                return False
        # Everything OK
        return True

    def nrp_maker(self,
        homology,
        seq_constr,
        struct_constr,
        target_size,
        background=None,
        part_type='RNA',
        struct_type=None,
        synth_opt=True,
        local_model_fn=None,
        global_model_fn=None,
        jump_count=10,
        fail_count=1000,
        output_file=None,
        verbose=True,
        abortion=True,
        allow_internal_repeat=False):

        if verbose:
            print '\n[Non-Repetitive Parts Calculator - Maker Mode]'
        build_parts = True

        # Check Maker Constraints
        if verbose:
            print '\n[Checking Constraints]'
            print '  Sequence Constraint: {}'.format(seq_constr)
            print ' Structure Constraint: {}'.format(struct_constr)
            print '    Target Size      : {} parts'.format(target_size)
            print '           Lmax      : {} bp'.format(homology-1)
            print '  Internal Repeats   : {}'.format(allow_internal_repeat)
        check_status = self._check_maker_constraints(
            seq_constr,
            struct_constr,
            allow_internal_repeat,
            target_size,
            homology)

        if check_status == False:
            if verbose:
                print ' Check Status: FAIL\n'
            build_parts = False
        else:
            if verbose:
                print '\n Check Status: PASS'
        
        # Background Check
        if build_parts:
            if not background is None:
                if verbose:
                    print '\n[Checking Background]\n Background: {}'.format(background)
                if isinstance(background, kmerSetDB):
                    self.background = background
                    if background.K != homology:
                        build_parts = False
                        print '\n [ERROR]    Background Lmax is {}, but Constraint Lmax is {}'.format(
                            background.K-1,
                            homology-1)
                        print ' [SOLUTION] Try correcting Lmax\n'
                        if verbose:
                            print ' Check Status: FAIL\n'
                    else:
                        if verbose:
                            print '\n Check Status: PASS'

                else:
                    build_parts = False
                    print '\n [ERROR]    Background Object is INVALID'
                    print ' [SOLUTION] Try instantiating background via nrpcalc.background(...)\n'
                    if verbose:
                        print ' Check Status : FAIL\n'

        # Arguments Check
        if build_parts:
            if verbose:
                print '\n[Checking Arguments]'
                print '   Part Type : {}'.format(part_type)
                print ' Struct Type : {}'.format(struct_type)
                print '  Synth Opt  : {}'.format(synth_opt)
                print '   Jump Count: {}'.format(jump_count)
                print '   Fail Count: {}'.format(fail_count)
                print ' Output File : {}'.format(output_file)
            check_status = self._check_maker_inputs(
                part_type,
                struct_type,
                synth_opt,
                jump_count,
                fail_count,
                output_file,
                verbose)

            if check_status == False:
                if verbose:
                    print ' Check Status: FAIL\n'
                build_parts = False
            else:
                if verbose:
                    print '\n Check Status: PASS'

        if not build_parts:
            # Cleanups
            self.background = None
            self.kmer_db.drop()
            self.kmer_db = None
            raise RuntimeError('Invalid Constraints, Background or Arguments')
        print

        # kmerSetDB Setup
        projector.setup_proj_dir(self.proj_id)
        self.kmer_db = kmerSetDB(
            path='./{}/kmerDB'.format(self.proj_id),
            homology=homology,
            verbose=verbose)

        # Project Setup
        if output_file is None:
            projector.setup_proj_dir(self.proj_id)
            output_file = './{}/seq_list.fa'.format(self.proj_id)

        # Execute Maker
        with open(output_file, 'w') as out_file:
            seq_constr    = seq_constr.upper().replace('U', 'T')
            struct_constr = self._get_adjusted_struct(
                struct_constr,
                seq_constr)
            current_nrp_count = 0
            memory_exhausted  = False

            if verbose:
                print 'Constructing Toolbox:\n'

            for non_coding_nrp in self._get_non_coding_nrps(
                homology,
                seq_constr,
                struct_constr,
                struct_type,
                target_size,
                synth_opt,
                local_model_fn,
                global_model_fn,
                jump_count,
                fail_count,
                verbose,
                abortion,
                allow_internal_repeat):

                # Write out genetic part
                if non_coding_nrp:
                    out_file.write(
                        '>non-repetitive part {}\n'.format(
                            current_nrp_count+1))
                    non_coding_nrp = '\n'.join(
                        textwrap.wrap(
                            non_coding_nrp, 80))
                    out_file.write('{}\n'.format(
                        non_coding_nrp))
                    current_nrp_count += 1
                else:
                    if non_coding_nrp is None:
                        print 'Failure Limits Exceeded or k-mers Exhausted. Cannot Build More Parts.'
                    else:
                        print 'Memory Capacity at Full. Cannot Build More Parts.'
                        memory_exhausted = True

                # Memory no longer available ... stop
                if memory_exhausted:
                    break
            if verbose:
                print '\nConstruction Complete.\n'

        # Detach Background
        self.background = None

        # Remove kmerSetDB
        self.kmer_db.drop()
        self.kmer_db = None

        # Pack output in dictionary
        parts_dict = {}
        for i,line in enumerate(utils.stream_fasta_seq_list(output_file)):
            line = line.strip()
            parts_dict[i] = line
       
        # Cleanups and Return
        projector.remove_proj_dir()
        if verbose:
            print 'Non-Repetitive Toolbox Size: {}'.format(current_nrp_count)
        return parts_dict

def main():
    homology = 16
    sm_obj = NRPMaker(seed=7)

    # Hammerhead Nielsen Paper
    seq         = 'NNNN AGNNNU CANNNNN UGUGCUU NNNNNU CUGAUGA NNNN GUGA NNNN GAAA NNNC CUCU NNNNN UAAU NNNNN UUAA NNNN' # Nielsen Like
    struct      = 'xxxx x((((( x(((((( xxxxxxx )))))) xxxxxxx (((( xxxx )))) xxx) )))) xxxx ((((( xxxx ))))) xxxx xxxx'
    seq         = ''.join(seq.split(' '))
    struct      = ''.join(struct.split(' '))
    output_file = None #'riboz.fa'
    background  = None #utils.get_fasta_seq_list(fasta_filename='input.fa.bk104')

    # Final Result Store
    final_toolbox = {}

    # Initialize Background
    background = kmerSetDB(
        path='./testDB',
        homology=homology,
        verbose=True)
    # background.multiadd(
    #     utils.get_fasta_seq_list(
    #         fasta_filename='input.fa.bk104'))

    # Background Based Single Part Design Works
    t0 = time()
    tt = 0
    toolbox1 = sm_obj.nrp_maker(homology, [seq], [struct], [10],
        struct_type='mfe',
        background=background,
        jump_count=100,
        fail_count=1000,
        synth_opt=False,
        verbose=True,
        abortion=True,
        allow_internal_repeat=True,
        output_file=output_file)
    tt += time() - t0
    final_toolbox.update(
        zip(range(len(final_toolbox), len(final_toolbox)+len(toolbox1)),
            toolbox1.values()))

    # Adding More Background Post Design Works
    background.multiadd(toolbox1.values())

    # Serial Designs from Multiple Constraints Works
    t0 = time()
    toolbox2 = sm_obj.nrp_maker(homology, [seq]*2, [struct]*2, [20]*2,
        struct_type='mfe',
        background=background,
        jump_count=100,
        fail_count=1000,
        synth_opt=False,
        verbose=True,
        abortion=True,
        allow_internal_repeat=False,
        output_file=output_file)
    tt += time() - t0
    final_toolbox.update(
        zip(range(len(final_toolbox), len(final_toolbox)+len(toolbox2)),
            toolbox2.values()))

    # Assert All Parts Unique and Non-Repetitive
    assert len(set(final_toolbox.values())) == len(
        finder.nrp_finder(final_toolbox.values(), homology, None, verbose=False))

    # Drop Background
    background.drop()

    # Report Time Elapsed
    print '\nWall Time {} sec'.format(tt)

if __name__ == '__main__':
    main()
