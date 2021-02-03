class kmerSetArray(object):
    '''
    Custom set-array structure for
    k-mer storage.
    '''

    def __init__(self, size):
        '''
        kmerSetArray constructor.
        '''
        try:
            self.unique_kmers = dict()
            self.indexd_kmers = [' '] * size
        except Exception as E:
            raise E

    def __setitem__(self, index, kmer):
        '''
        User function to set kmer at given index.
        '''
        try:
            if not self.indexd_kmers[index] == ' ':
                if self.indexd_kmers[index] in self.unique_kmers:
                    self.unique_kmers.pop(
                        self.indexd_kmers[index])
            self.unique_kmers[kmer]  = index
            self.indexd_kmers[index] = kmer
        except Exception as E:
            raise E

    def __contains__(self, kmer):
        '''
        User function to check existence of kmer
        in kmerSetArray.
        '''
        if kmer in self.unique_kmers:
            return True
        else:
            return False

    def index(self, kmer):
        '''
        User function to return the index of
        kmer if stored in kmerSetArray.
        '''
        if kmer in self.unique_kmers:
            return self.unique_kmers[kmer]
        return None

    def __len__(self):
        '''
        User function to return the number of
        kmers stored in kmerSetArray.
        '''
        return len(self.unique_kmers)

    def __repr__(self):
        '''
        User function to return a string
        representation of kmerSetArray.
        '''
        return 'kmerSetArray[{}]'.format(
            ', '.join(map(str, self.indexd_kmers)))

    def __str__(self):
        '''
        See __repr__
        '''
        return self.__repr__()

    def __del__(self):
        '''
        kmerSetArray destructor.
        '''
        self.unique_kmers = None
        self.indexd_kmers = None

def test():

    # required test imports
    import random
    import time

    from . import utils

    # returns a random DNA string of length
    def get_DNA(length):
        return ''.join(random.choice('ATGC') for _ in xrange(length))

    # create seq
    seq = get_DNA(length=100)
    print('Testing Operations on {}...'.format(seq[:25], seq[-25:]))

    # object test
    mySet = kmerSetArray(size=len(seq)-16+1)

    # __setitem__ test
    t0 = time.time()
    for i,kmer in enumerate(utils.stream_kmers(seq, k=16)):
        mySet[i] = kmer
    print('__setitem__  time = {} sec'.format(time.time()-t0))

    # __contains__ test
    t0 = time.time()
    for i,kmer in enumerate(utils.stream_kmers(seq, k=16)):
        assert kmer in mySet
    print('__contains__ time = {} sec'.format(time.time()-t0))

    # index test
    t0 = time.time()
    for i,kmer in enumerate(utils.stream_kmers(seq, k=16)):
        assert mySet.index(kmer) == i
    print('index        time = {} sec'.format(time.time()-t0))

    # __len__ test
    t0 = time.time()
    assert len(mySet) == len(seq)-16+1
    print('__len__      time = {} sec'.format(time.time()-t0))

    # __del__ test
    t0 = time.time()
    del mySet
    print('__del__      time = {} sec'.format(time.time()-t0))


if __name__ == '__main__':
    test()