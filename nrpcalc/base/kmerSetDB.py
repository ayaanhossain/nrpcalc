import sys
import shutil

import utils
import plyvel


class kmerSetDB(object):
    '''
    Plyvel based scalable on disk kmerSetDB
    for k-mer storage.
    '''
    
    def __init__(self, path, homology, verbose=False):
        '''
        kmerSetDB constructor.
        '''
        try:
            # Assert path is string
            if not isinstance(path, str):
                print '\n[Non-Repetitive Parts Calculator - Background]'
                print '\n [ERROR]    Path must be a string, not {}'.format(
                    path)
                print ' [SOLUTION] Try correcting Path\n'
                raise ValueError

            # Assert homology is positive integer
            if not isinstance(homology, int):
                print '\n[Non-Repetitive Parts Calculator - Background]'
                print '\n [ERROR]    Lmax must be an integer, not \'{}\''.format(homology-1)
                print ' [SOLUTION] Try correcting Lmax\n'
                raise ValueError
            if homology-1 < 5:
                print '\n[Non-Repetitive Parts Calculator - Background]'
                print '\n [ERROR]    Lmax must be greater than 4, not \'{}\''.format(homology-1)
                print ' [SOLUTION] Try correcting Lmax\n'
                raise ValueError
            
            # Path setup
            self.PATH = path.rstrip('/') + '/'

            # Create/Open Plyvel object
            self.DB = plyvel.DB(
                name=self.PATH,
                create_if_missing=True,
                error_if_exists=False)

            # Length setup
            if self.DB.get('LEN') is None:
                self.DB.put('LEN', '0')
                self.LEN = 0
            else:
                self.LEN = int(self.DB.get('LEN'))

            # K setup
            if self.DB.get('K') is None:
                self.DB.put('K', str(homology))
                self.K = int(homology)
            else:
                self.K = int(self.DB.get('K'))

            # Verbosity setup
            self.VERB = bool(verbose)
            
            # Object ALIVE status
            self.ALIVE = True

        except ValueError as E:
            raise ValueError('Invalid Path or Lmax')
        except Exception as E:
            raise E

    def __repr__(self):
        '''
        User function to return the string representation of
        kmerSetDB.
        '''
        return 'kmerSetDB stored at {} with {} {}-mers'.format(
            self.PATH, self.LEN, self.K)

    def __str__(self):
        '''
        See __repr__
        '''
        return self.__repr__()

    def alivemethod(method):
        '''
        Internal decorator to gate kmerSetDB operation
        once dropped.
        '''
        def wrapper(self, *args, **kwargs):
            if self.ALIVE:
                return method(self, *args, **kwargs)
            else:
                raise RuntimeError('kmerSetDB was closed or dropped')
        return wrapper

    def __len__(self):
        '''
        User function to return the number of keys stored
        in kmerSetDB.
        '''
        return self.LEN

    def _verb_action(self, action, index, seq):
        '''
        Internal helper function to print sequence
        updates to kmerSetDB based on verbosity.
        '''
        if self.VERB:
            psl = min(len(seq), 10)
            printed = '  {} Seq {}: {}...'.format(
                action, index, seq[:psl])
            clear_length = len(printed)
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            sys.stdout.flush()

    @alivemethod
    def add(self, seq):
        '''
        User function to add seq k-mers to kmerSetDB.
        '''
        if not seq in ['K', 'LEN']:
            try:
                with self.DB.write_batch(
                    transaction=True,
                    sync=False) as wb:
                    for kmer in utils.stream_min_kmers(
                        seq=seq, k=self.K):
                        if self.DB.get(kmer) is None:
                            wb.put(kmer, '1')
                            self.LEN += 1
                    self.DB.put('LEN', str(self.LEN))
            except Exception as E:
                raise E

    @alivemethod
    def multiadd(self, seq_list):
        '''
        User function to add seq k-mers for each seq in
        seq_list to kmerSetDB via single transaction.
        '''
        try:
            with self.DB.write_batch(
                transaction=True,
                sync=False) as wb:
                pt = False
                for index, seq in enumerate(seq_list):
                    if not pt and self.VERB:
                        print '\n[Background Processing]'
                        pt = True
                    if not seq in ['K', 'LEN']:
                        self._verb_action(
                            action='Adding',
                            index=index,
                            seq=seq)
                        for kmer in utils.stream_min_kmers(
                            seq=seq, k=self.K):
                            if self.DB.get(kmer) is None:
                                wb.put(kmer, '1')
                                self.LEN += 1
                self.DB.put('LEN', str(self.LEN))
                if self.VERB:
                    print
        except Exception as E:
            raise E

    @alivemethod
    def __contains__(self, seq):
        '''
        Python dunder function to check existence of
        any k-mer from seq in kmerSetDB.
        '''
        try:
            for kmer in utils.stream_min_kmers(
                seq=seq, k=self.K):
                if self.DB.get(kmer) == '1':
                    return True
            return False
        except Exception as E:
            raise E

    def multicheck(self, seq_list):
        '''
        User function to check existence of any k-mer for
        each seq in seq_list in kmerSetDB via single
        transaction.
        '''
        try:
            for seq in seq_list:
                for kmer in utils.stream_min_kmers(
                    seq=seq, k=self.K):
                    if self.DB.get(kmer) == '1':
                        yield True
                        break
                else:
                    yield False
        except Exception as E:
            raise E

    def __iter__(self):
        '''
        User fuction to iterate over k-mers stored in
        kmerSetDB.
        '''
        with self.DB.iterator() as it:
            for key, _ in it:
                if not key in ['K', 'LEN']:
                    yield key

    def __len__(self):
        '''
        User function to check the number of k-mers stored
        in kmerSetDB.
        '''
        return self.LEN

    @alivemethod
    def remove(self, seq):
        '''
        User function to remove all k-mers in seq
        from kmerSetDB.
        '''
        if not seq in ['K', 'LEN']:
            try:
                with self.DB.write_batch(
                    transaction=True,
                    sync=False) as wb:
                    for kmer in utils.stream_min_kmers(
                        seq=seq, k=self.K):
                        if self.DB.get(kmer) == '1':
                            wb.delete(kmer)
                            self.LEN -= 1
                    self.DB.put('LEN', str(self.LEN))
            except Exception as E:
                raise E

    @alivemethod
    def multiremove(self, seq_list):
        '''
        User function to remove all k-mers from each
        seq in seq_list from kmerSetDB via single
        transaction.
        '''
        try:
            with self.DB.write_batch(
                transaction=True,
                sync=False) as wb:
                pt = False
                for index, seq in enumerate(seq_list):
                    if not pt and self.VERB:
                        print '\n[Background Processing]'
                        pt = True
                    if not seq in ['K', 'LEN']:
                        self._verb_action(
                            action='Removing',
                            index=index,
                            seq=seq)
                        for kmer in utils.stream_min_kmers(
                            seq=seq, k=self.K):
                            if self.DB.get(kmer) == '1':
                                wb.delete(kmer)
                                self.LEN -= 1
                self.DB.put('LEN', str(self.LEN))
                if self.VERB: print
        except Exception as E:
            raise E

    def clear(self):
        '''
        User function to clear all k-mers stored in
        kmerSetDB.
        '''
        if self.LEN:
            self.multiremove(iter(self))

    def close(self):
        '''
        User function to close kmerSetDB.
        '''
        if self.ALIVE:
            self.DB.close()
            self.DB = None
            self.ALIVE = False
            return True
        return False

    def drop(self):
        '''
        User function to drop kmerSetDB.
        '''
        if self.ALIVE:
            self.DB.close()
            self.DB = None
            shutil.rmtree(self.PATH)
            self.ALIVE = False
            return True
        return False

def test():

    # required test imports
    import random
    import time

    # returns a random DNA string of length
    def get_DNA(length):
        return ''.join(random.choice('ATGC') for _ in xrange(length))
    
    # create corpus
    total_elements = 200000
    dna_strings = sorted(get_DNA(100) for _ in xrange(total_elements))
    print 'Testing Operations on 400000 items'

    # clearance
    try:
        shutil.rmtree('./testBG/')
    except:
        pass

    # object test
    myBG = kmerSetDB(homology=99, path='./testBG')
    assert isinstance(myBG, kmerSetDB)
    assert repr(myBG) == 'kmerSetDB stored at ./testBG/ with 0 99-mers'

    # add test
    t0 = time.time()
    for i in xrange(total_elements):
        myBG.add(dna_strings[i])
    print 'add          time = {} sec'.format(time.time()-t0)
    myBG.clear()

    # multiadd test
    t0 = time.time()
    myBG.multiadd(dna_strings)
    print 'multiadd     time = {} sec'.format(time.time()-t0)

    # __contains__ test
    assert 'Uninserted String' not in myBG
    t0 = time.time()
    for dna in dna_strings:
        assert dna in myBG
    print '__contains__ time = {} sec'.format(time.time()-t0)

    # multicheck test
    t0 = time.time()
    assert all(myBG.multicheck(seq_list=dna_strings))
    print 'multicheck   time = {} sec'.format(time.time()-t0)

    # __iter__ test
    t0 = time.time()
    sorted(myBG)
    print '__iter__     time = {} sec'.format(time.time()-t0)

    # __len__ test
    t0 = time.time()
    assert len(myBG) == total_elements * 2
    print '__len__      time = {} sec'.format(time.time()-t0)

    # remove test
    t0 = time.time()
    for i in xrange(total_elements):
        myBG.remove(dna_strings[i])
    print 'remove       time = {} sec'.format(time.time()-t0)
    assert len(myBG) == 0

    # multiremove test
    myBG.multiadd(dna_strings)
    t0 = time.time()
    myBG.multiremove(dna_strings)
    print 'multiremove  time = {} sec'.format(time.time()-t0)
    assert len(myBG) == 0

    # clear test
    myBG.multiadd(dna_strings)
    t0 = time.time()
    myBG.clear()
    print 'clear        time = {} sec'.format(time.time()-t0)

    # close test
    t0 = time.time()
    myBG.close()
    print 'close        time = {} sec'.format(time.time()-t0)

    # drop test
    myBG = kmerSetDB(homology=99, path='./testBG')
    t0 = time.time()
    myBG.drop()
    print 'drop         time = {} sec'.format(time.time()-t0)

if __name__ == '__main__':
    test()