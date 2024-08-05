class ShuffleAlgebraImproved(sage.algebras.shuffle_algebra.ShuffleAlgebra):
    """
    This class is a reimplementation of Sage's built-in ShuffleAlgebra class, with two improvements:
    1) 3x efficiency saving in the built-in to_dual_pbw_element function.
    2) Enable parallelization of the to_dual_pbw_element function.
    """
    
    def __init__(self, R, names, prefix=None):
        # Will be set to true if to_dual_pbw_element has been calculated for all relevant words 
        self.precomputed = False 
        super().__init__(R, names, prefix)

    @staticmethod
    def __classcall_private__(cls, R, names, prefix=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = ShuffleAlgebra(QQ, 'xyz')
            sage: F2 = ShuffleAlgebra(QQ, ['x','y','z'])
            sage: F3 = ShuffleAlgebra(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        if prefix is None:
            prefix = 'B'
        return super().__classcall__(cls, R,
                                     Alphabet(names), prefix)

    def to_dual_pbw_element(self, w):
        """
            A rewriting of the built in function ShuffleAlgebra.to_dual_pbw_element.
            A simple change leads to a 3x speedup
        """
        if self.precomputed:
            return self.cache[w]
        
        D = self.dual_pbw_basis()
        l = {}
        
        while w != self.zero():
            (min_elt, coeff) = max(w, key=lambda x: x[0])
            l[min_elt] = l.get(min_elt, 0) + coeff
            w -= coeff * D.expansion_on_basis(min_elt)
        
        return D.sum_of_terms((m, c) for m, c in l.items() if c != 0)

    def to_dual_pbw_element_async(self, w, cache):
        """
            A rewriting of the built-in function ShuffleAlgebra.to_dual_pbw_element 
            to enable parallelization. The only difference is that the cache
            that is normally used by expansion on basis is now handled manually
            rather than using the @cached_method decorator
        """
        original_w = w
        D = self.dual_pbw_basis()
        l = {}
        
        while w != self.zero():
            (min_elt, coeff) = max(w, key=lambda x: x[0])
            l[min_elt] = l.get(min_elt, 0) + coeff
            w -= coeff * self.expansion_on_basis_async(min_elt, cache)
        
        return original_w, D.sum_of_terms((m, c) for m, c in l.items() if c != 0)

    def expansion_on_basis_async(self, w, cache):
        """
            A copy of the built-in expansion_on_basis function to enable parallelization
        """
        if w in cache:
            return cache[w]
    
        if not w:
            return self.one()
        if len(w) == 1:
            return self.monomial(w)
        if w.is_lyndon():
            W = self.dual_pbw_basis().basis().keys()
            letter = W([w[0]])
            expansion = self.expansion_on_basis_async(W(w[1:]),cache)
    
            cache[w] = self.sum_of_terms((letter * i, c)
                                          for i, c in expansion)
    
            return cache[w]
            
        lf = w.lyndon_factorization()
        powers = {}
        for i in lf:
            powers[i] = powers.get(i, 0) + 1
        denom = prod(factorial(p) for p in powers.values())
        result = self.prod(self.expansion_on_basis_async(i,cache) for i in lf)
        cache[w] = self(result / denom)
        return cache[w]
