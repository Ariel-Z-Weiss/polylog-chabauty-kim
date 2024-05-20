import string
from functools import cache
from itertools import groupby
import numpy
import cProfile
import pstats
from multiprocessing import Pool, Manager


class ShuffleAlgebraPolynomialRing(sage.symbolic.ring.SymbolicRing):
    """
    Constructor for the ShuffleAlgebraPolynomialRing class, a subclass of sage.symbolic.ring.SymbolicRing.
    
    This was originally a subclass of MPolynomialRing_libsingular, however, it seems the polynomial rings are limited to 2^15
    variables, which we will need to exceed.

    Let S be a finite set of primes. The ring O(U_S) is a shuffle algebra on words in the alphabet
    - \tau_p for p in S
    - \sigma_{2n + 1} for integers n \ge 1.

    If we assign each \tau_p degree 1, and each \sigma_{2n+1} degree 2n+1, then O(U_S) is a graded ring.

    As a shuffle algebra over Q, O(U_S) is isomorphic to a polynomial ring, whose generators are Lyndon 
    words. Let d = 'halfweight'. Then, by Definition 2.6, the ring O(U_S)_{\le d} is the subring generated by words
    degree at most d. In practice, we only need to consider those words in which each \sigma_{2n+1} appears
    at most once.

    An object of this class is a symbolic ring over Q whose generators are the Lyndon words of degree
    at most d, and such that each \sigma_{2n+1} appears at most once.

    Parameters:
    halfweight (int): The highest degree of a Lyndon word generator of the ring.
    number_of_primes (list): The number of prime numbers in the set S.
"""

    def __init__(self, halfweight, number_of_primes):
        self.halfweight = halfweight
        self.number_of_primes = number_of_primes
        
        self.sigma_letters = self._sigma_letters_and_degrees(halfweight)
        self.tau_letters = self._tau_letters_and_degrees(number_of_primes)
        self.letters = [letter[1] for letter in self.sigma_letters] + [letter[1] for letter in self.tau_letters]
        
        self.variables = self._generate_lyndon_words(
            self.halfweight,
            self.sigma_letters,
            [letter[1] for letter in self.tau_letters]
        )
        
        self.number_of_variables = len(self.variables)
        
        self.Words = Words(self.letters)

        # Initialise the parent class as a ring over Q
        super().__init__(QQ)    

        # Add all our variables so that they can be symbolically manipulated
        self.var(self.variables)
        
    def _sigma_letters_and_degrees(self, halfweight):
        """
        Return the list of letters \sigma_{3}, \sigma_{5}, ... up to \sigma_{halfweight} or \sigma_{halfweight - 1}
        depending on the parity of 'halfweight' along with their associated degrees.
        
        A bug in the sage shuffle algebra package means that 'sigma_3' would be encoded as a workd of length
        7. To bypass this bug, we encode \sigma_{2n+1} by the n-th capital letter.
        """
        return [(f'sigma{2*n + 1}', 
                 string.ascii_uppercase[n-1],
                 2*n + 1) 
                for n in range(1, floor((halfweight - 1)/2) + 1)] 
    
    def _tau_letters_and_degrees(self, number_of_primes):
        """
        Return the list of letters \tau_p for each prime p in S along with their associated degrees.
        
        A bug in the sage shuffle algebra package means that 'tau_p' would be encoded as a workd of length
        7. To bypass this bug, if we encode \tau_{p_n} as the n-th lowercase letter.
        """
        return [(f'tau{index}',string.ascii_lowercase[index], 1)
                for index in range(number_of_primes)]
    
    @staticmethod
    def words_of_fixed_length(length, letters):
        """
        Takes as input an ordered list of letters, and outputs the words in those letters of given length.
        
        Inputs:
        - length: the length of each word
        - tau_letters: a list of letters
                       
        Output:
        - a list of strings consisting of the words
        """
        for word in Words(letters, length):
            yield str(word)
    
    @staticmethod
    def lyndon_words_of_fixed_length(length, letters):
        """
        Takes as input an ordered list of letters, and outputs the Lyndon words in those letters and
        given length.

        This code was taken from https://www.geeksforgeeks.org/python-program-for-generating-lyndon-words-of-length-n/
        
        Inputs:
        - length: the length of each Lyndon word
        - tau_letters: a list of letters
                       
        Output:
        - a list of strings consisting of the Lyndon words
        """
        number_of_letters = len(letters)
        # To store the indices 
        # of the characters 
        w = [-1] 

        # Loop till w is not empty 
        while w: 

            w[-1] += 1
            m = len(w) 
            if m == length: 
                yield 'S' + ''.join(letters[i] for i in w)
                
            # Repeating w to get a 
            # n-length string 
            while len(w) < length: 
                w.append(w[-m]) 
            
            # Removing the last character 
            # as long it is equal to 
            # the largest character in S 
            while w and w[-1] == number_of_letters - 1: 
                w.pop()
            
        return
     
    def _generate_lyndon_words_of_fixed_length(self, length, sigma_letters_and_degrees, tau_letters):
        """
        Return a list of all the Lyndon words in our alphabet, of given, and containing
        at most one sigma letter.
        
        Inputs:
        - length (int): the length  
        - sigma_letters_and_degrees: a list of tuples (letter_string, encoded_letter_string, degree)
        as produced by self._sigma_letters_and_degrees
        - tau_letters: a list of letters
        """
        
        # A Lyndon word of length n is either
        # - a Lyndon word of length n in the tau_letters
        # - \sigma_{2k+1} followed by *any* word of length n - (2k+1) in the tau letters
      
        # Add the Lyndon words that include a sigma
        for _, letter, degree in sigma_letters_and_degrees:
            remaining_length = length - degree
            if remaining_length < 0:
                break
                
            for word in ShuffleAlgebraPolynomialRing.words_of_fixed_length(remaining_length, tau_letters):
                yield "".join(['S', letter, word])
             
        # Add the Lyndon words not including a sigma
        yield from ShuffleAlgebraPolynomialRing.lyndon_words_of_fixed_length(length,tau_letters)
            
        
    def _generate_lyndon_words(self, max_length, sigma_letters_and_degrees, tau_letters):
        """
        Return a list of all the Lyndon words in our alphabet, of length up to 'max_length', and containing
        at most one sigma letter.
        
        Inputs:
        - max_length (int): the maximum length  
        - sigma_letters_and_degrees: a list of tuples (letter_string, encoded_letter_string, degree)
        as produced by self._sigma_letters_and_degrees
        - tau_letters: a list of letters
        """
        lyndon_words = []
        
        for length in range(1, max_length + 1):
            number_of_words = len(lyndon_words)
            lyndon_words.extend(self._generate_lyndon_words_of_fixed_length(
                            length,
                            sigma_letters_and_degrees,
                            tau_letters))

        return lyndon_words
    
    @cache
    def decompose_pbw_word(self, word):
        """
        Takes as input a word in the dual_pbw_basis of the shuffle algebra. 
        Decomposes that word as a product of S_w's where w is a Lyndon word.
        
        The output is a tuple of word strings
        """
        #Convert word from a ShuffleAlgebra object to a word object
        word = Word(self.Words(str(word)))
        return word.lyndon_factorization()
    
    def _convert_word_to_polynomial(self ,word):
        '''
        Takes 'word', a word in a ShuffleAlgebra, and converts it to a polynomial
        in this ring
        '''
        decomposed_word = self.decompose_pbw_word(word)
        polynomial = 1

        for lyndon_word, group in groupby(decomposed_word):
            exponent = len(list(group))
            polynomial *= self(f'S{str(lyndon_word)}')^exponent / factorial(exponent)
        return polynomial
    
    def convert_from_shuffle_algebra(self, sentence):
        '''
        Takes 'sentence', an element of a ShuffleAlgebra object, and converts it to a monomial
        in this ring
        '''
        polynomial = 0
        
        for word, coefficient in sentence:
            polynomial += coefficient * self._convert_word_to_polynomial(word)
        
        return polynomial

### For some reason, python won't regonise sage.algebras.shuffle_algebra.ShuffleAlgebra
### until I call this. And for other strange reasons, I get an attribute error if I try to 
### just inherit from ShuffleAlgebra
no_idea_why = ShuffleAlgebra(QQ, 'ab')

class ShuffleAlgebraImproved(sage.algebras.shuffle_algebra.ShuffleAlgebra):
    """
    This class is a reimplemantation of Sage's built in ShuffleAlgebra class, in order to include two
    improvements:
    1) a major efficiency saving in the built in to_dual_pbw_element function
    2) To enable parallelisation of the to_dual_pbw_element function.
    """
    
    def __init__(self, R, names, prefix=None):
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
            A rewriting of the built in function ShuffleAlgebra.to_dual_pbw_element 
            to enable parallelisation. The only difference is that the cache
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
            A copy of the built in expansion_on_basis function to enable parallelisation
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
    
class PhiPolynomialRing(sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular):
    """
    Constructor for the PolylogPolynomialRing class, a subclass of MPolynomialRing_libsingular.
        
    An object of this class is a polynomial ring 
        Q[\Phi],
    as defined in Section 4.1. 
    
    Parameters:
    halfweight (int): The largest Li_i that we want to work with.
    number_of_primes (list): The number of prime numbers in the set S.
    
    """
    def __init__(self, halfweight, number_of_primes):
        self.halfweight = halfweight
        self.number_of_primes = number_of_primes
        
        self.variables, self.degrees = self._generate_phi_variables_and_degrees(halfweight)
        
        super().__init__(
            QQ,
            len(self.variables),
            self.variables,
            order=TermOrder('wdegrevlex', self.degrees)
        )
        
            
    def _generate_phi_variables_and_degrees(self, maximum_halfweight):
        '''
        Compute the phi variables \Phi_{\lambda}^w along with their degrees.
        
        The degrees must be alphanumeric.
        - we use phi0t{index} and phi1t{index} to represent \Phi_{e_0}^{\tau_p}
        - we use phisigma{i} to representation \Phi_{e_0e_1...e_1}^{\sigma_i}
        
        Returns a list of variables and a list of degrees
        '''
        phi_variables = []
        degrees = []

        if maximum_halfweight <= 0: return [], []
        
        # Add the variables \Phi_{e_0}^{\tau_p} and \Phi_{e_1}^{\tau_p} for p in S
        
        for i in range(2):
            for index in range(self.number_of_primes):
                phi_variables.append(f'phi{i}t{index}')
                degrees.append(1)
        
        for halfweight in range(2, maximum_halfweight + 1):
            #There are no new phi variables in odd halfweights
            if halfweight % 2 == 0:
                continue
                
            phi_variables.append(f'phisigma{halfweight}')
            degrees.append(halfweight)
            
        return phi_variables, degrees
    
    def monomials_of_degree(self, degree):
        '''
        Return a list of all monomials of given degree. 
        '''
        weighted_vectors = list(WeightedIntegerVectors(degree, self.degrees))
        
        return [prod([self(variable)^degrees[index] for index, variable in enumerate(self.variables)])
               for degrees in weighted_vectors]

class ThetaSharpOperator:
    """
    An object of this class represents the operator
        Theta^sharp: \O(\Pi_d) \otimes O(U_S)_{\le d} \to O(U_S)_{\le d}[\Phi, d]
    acting on monomials in 'log' and 'Li_i'.
    """
    
    def __init__(self, halfweight, number_of_primes):
        self.halfweight = halfweight
        self.number_of_primes = number_of_primes

        #Create the shuffle algebra O(U_S)_{\le d}
        self.OU_S_ring = ShuffleAlgebraPolynomialRing(halfweight, number_of_primes)

        self.tau_letters = [letter[1] 
                       for letter in self.OU_S_ring.tau_letters]
        self.sigma_letters = [letter[1] 
                       for letter in self.OU_S_ring.sigma_letters]
        self.letters = self.sigma_letters + self.tau_letters
        self.shuffle_algebra = ShuffleAlgebraImproved(QQ, self.letters)
        
        #Create the polynomial ring Q[\Phi, d]. 
        self.phi_algebra = PhiPolynomialRing(halfweight, number_of_primes)

        #Create the polynomial ring O(U_S)_{\le d}[\Phi, d]
        self.OU_phi_algebra = PolynomialRing(self.OU_S_ring, 
                                             len(self.phi_algebra.variables), 
                                             self.phi_algebra.variables,
                                             order=TermOrder('wdegrevlex', self.phi_algebra.degrees))
        
        # Create a dictionary mapping the letter corresponding to each \tau_p and \sigma_{2n+1}
        # to the corresponding generator \phi_{e_0}^\tau_p or \phi_{e_0....e1}^\sigma
        self.tau_letter_to_phi_dict = {letter : self.phi_algebra(f'phi0t{index}') 
                              for index,letter in enumerate(self.tau_letters)}
        self.sigma_letter_to_phi_dict = {letter : self.phi_algebra(f'phisigma{2*index+3}') 
                              for index,letter in enumerate(self.sigma_letters)}
        self.tau_letter_to_phi_e1_dict = {letter : self.phi_algebra(f'phi1t{index}') 
                              for index,letter in enumerate(self.tau_letters)}
        
        self.li_values_dict = {}
        self.li_values_random_dict = {}

    def generate_words_to_precompute(self, halfweight):
        """
            In the parallel version, we first precompute to_dual_pbw_element on every
            word that can appear in some Li(n) up to the halfweight
        """
        for n in range(1, halfweight + 1):
            for tau_letter in self.tau_letters:
                for word in ShuffleAlgebraPolynomialRing.words_of_fixed_length(n-1, self.tau_letters):
                    yield self.shuffle_algebra(tau_letter + word)
                    
            for sigma_letter_and_degree in self.OU_S_ring.sigma_letters:
                remaining_degree = n - sigma_letter_and_degree[2]
                if remaining_degree < 0:
                    continue
                
                sigma_letter = sigma_letter_and_degree[1]
                
                for word in ShuffleAlgebraPolynomialRing.words_of_fixed_length(remaining_degree, self.tau_letters):
                    yield self.shuffle_algebra(sigma_letter + word)
               
    
    def set_random_integers(self, clear=True):
        '''
        To enable consistent test data, assign once and for all a set of integer values on which
        to evaluate the polynomials
        '''
        if clear:
            # Clear any saved random li values
            self.li_values_random_dict.clear()
        
        #random_ints = [i + 1 for i in range(self.OU_S_ring.number_of_variables)]
        random_ints = [(-1)^i for i in range(self.OU_S_ring.number_of_variables)]

        return {gen : gen.count('a') + 1 for gen in self.OU_S_ring.variables}
        # return {gen : random_int for gen, random_int in zip(self.OU_S_ring.gens(), random_ints)}

    def generate_random_integers(self, clear=True, integer_range=100):
        '''
        For each generator of sef.OU_S_ring, generate a random integer
        
        The point of including a gcd is to ensure that we are always working with integer arithmetic
        
        '''
        if clear:
            # Clear any saved random li values
            self.li_values_random_dict.clear()
        
        random_ints = numpy.random.randint(low=1, high=integer_range, size=self.OU_S_ring.number_of_variables)
        return {gen : Integer(random_int) for gen, random_int in zip(self.OU_S_ring.variables, random_ints)}
        
    @staticmethod
    def evaluate(polynomial, eval_dict):
        """
        Evaluates a polynomial at the entries of a dictionary, eval_dict.

        Previous iterations of this code used the built in .subs() method, which is very expensive!

        Returns: the evaluation of the polynomial, as a rational number.
        """
        polynomial = polynomial.polynomial(QQ)             
        evaluation = 0
        
        if len(polynomial.variables()) == 1:
            variable = polynomial.variables()[0]
            for exponent, coeff in enumerate(polynomial):
                evaluation += coeff * eval_dict[str(variable)]^exponent

            return evaluation
        
        for coefficient, monomial in polynomial:
            evaluation += coefficient * ThetaSharpOperator.evaluate_monomial(monomial, eval_dict)

        return evaluation

    @staticmethod
    def evaluate_monomial(monomial, eval_dict):
        """
        Evaluates a monomial at the entries of a dictionary, eval_dict.

        Returns: the evaluation of the polynomial, as a rational number.
        """
        evaluation = 1

        for variable in monomial.variables():
            evaluation *= eval_dict[str(variable)] ^ monomial.degree(variable)

        return evaluation
        
    def log(self, random_evaluation=False, eval_dict=None):
        """
        Returns the the image of log under theta_sharp. The output is a PolynomialDict 
        object whose coefficients are elements of a ShuffleAlgebraPolynomialRing object
        and whose monomials are elements of a PhiPolynomialRing object.
        
        If random_evaluation=True, compute the coefficient as the rational number obtained by
        randomly evaluating the corresponding shuffle algebra element
        
        Returns:
        PolynomialDict: The logarithm polynomial.
        """
        if random_evaluation and not eval_dict:
            raise ValueError('Evaluation dictionary not passed to function')

        if random_evaluation:
            algebra = self.phi_algebra

        else:
            algebra = self.OU_phi_algebra
       
        log_value = algebra(0)
        
        monomials_of_degree = algebra.monomials_of_degree(1)
        monomials_of_degree.reverse()
        # theta^sharp(log) is the sum of f_{\tau_p}\Phi_{e_0}^{\tau_p} for p in S
        for index, monomial in enumerate(monomials_of_degree):
            # Throw out the \Phi_{e_1}'s'
            if index >= self.number_of_primes: break
            
            #Pick out the \tau_p corresepnding to \Phi_{e_0}^{\tau_p}
            word = self.OU_S_ring.tau_letters[index][1]
            coefficient = self.OU_S_ring(f'S{word}')
            
            if random_evaluation:
                coefficient = ThetaSharpOperator.evaluate(coefficient, eval_dict)
                #coefficient = QQ(coefficient.subs(eval_dict))
                
            log_value += coefficient * monomial
            
        #Store the value
        if random_evaluation:
            self.li_values_random_dict['log'] = log_value

        else:
            self.li_values_dict['log'] = log_value
            
        return log_value
    
    def Li(self, n, random_evaluation=False, eval_dict=None):
        """
        Returns the the image of Li_n under theta_sharp. The output is a PolynomialDict 
        object whose coefficients are elements of a ShuffleAlgebraPolynomialRing object
        and whose monomials are elements of a PhiPolynomialRing object.
        
        If random_evaluation=True, compute the coefficient as the rational number obtained by
        randomly evaluating the corresponding shuffle algebra element.
        
        Returns:
        PolynomialDict: The Li(n) polynomial.
        """
        if random_evaluation and not eval_dict:
            raise ValueError('Evaluation dictionary not passed to function')


        #If the value has already been computed, do not repeat the computation
        if random_evaluation:
            if f'Li{n}' in self.li_values_random_dict:
                return self.li_values_random_dict[f'Li{n}']

            else:
                algebra = self.phi_algebra

        else:
            if f'Li{n}' in self.li_values_dict:
                return self.li_values_dict[f'Li{n}']

            else:
                algebra = self.OU_phi_algebra
       
        li_value = algebra(0)
        
        #Compute the coefficients that don't have any sigma variables
        #There is a coefficient for each pair consisting of a tau variable and a word in the
        #tau varibles of length n-1.
        for tau_letter in self.tau_letters:
            for word in ShuffleAlgebraPolynomialRing.words_of_fixed_length(n-1, self.tau_letters):
                #The coffiecient is f_{tau_letter + word}. Convert this word into the Lyndon basis, so it 
                #can be viewed as an element of self.OU_S_ring.
                sentence = self.shuffle_algebra.to_dual_pbw_element(self.shuffle_algebra(tau_letter + word))
                coefficient = self.OU_S_ring.convert_from_shuffle_algebra(sentence)

                if random_evaluation:
                    coefficient = ThetaSharpOperator.evaluate(coefficient, eval_dict)
                    #coefficient = QQ(coefficient.subs(eval_dict))
                #The variable is prod(\Phi_{e_0}^{\tau_p} : tau_p in word)\Phi_{e_1}^{tau_letter} 
                monomial = prod([algebra(self.tau_letter_to_phi_dict[letter]) for letter in word]) 
                monomial *= algebra(self.tau_letter_to_phi_e1_dict[tau_letter])

                li_value += coefficient * monomial
                
        #Compute the coefficients that have a sigma variable sigma_{2k+1}
        #There is a coefficient for each word in the tau varibles of length k - (2n+1).
        for sigma_letter_and_degree in self.OU_S_ring.sigma_letters:
            remaining_degree = n - sigma_letter_and_degree[2]
            if remaining_degree < 0:
                continue
            
            sigma_letter = sigma_letter_and_degree[1]
            
            for word in ShuffleAlgebraPolynomialRing.words_of_fixed_length(remaining_degree, self.tau_letters):
                #The coffiecient is f_{sigma_letter + word}. Convert this word into the Lyndon basis, so it 
                #can be viewed as an element of self.OU_S_ring.
                sentence = self.shuffle_algebra.to_dual_pbw_element(self.shuffle_algebra(sigma_letter + word))
                coefficient = self.OU_S_ring.convert_from_shuffle_algebra(sentence)
                
                if random_evaluation:
                    coefficient = ThetaSharpOperator.evaluate(coefficient, eval_dict)
                    #coefficient = QQ(coefficient.subs(eval_dict))
                
                #The variable is the product of \Phi_{e_0}^{\tau_p} for \tau_p appearing in the word,
                #multiplied by \Phi_{e_0e_1...e_1}^{\sigma}
                monomial = prod([algebra(self.tau_letter_to_phi_dict[letter]) for letter in word])
                monomial *= algebra(self.sigma_letter_to_phi_dict[sigma_letter])

                li_value += coefficient * monomial
            
        # Pythons .left_nullity() runs significantly faster on integer matrices than
        # on rational ones. On the other hand, multiplying each Li(n) by a scalar depending
        # on n will not affect the dimension of the kernel. So we can clear denominators
        # at this stage.
        if random_evaluation:
            values = [QQ(value) for value in li_value.coefficients()]
            gcd_of_values = gcd(values)

            li_value = li_value / gcd_of_values
            
        #Store the value
        if random_evaluation:
            self.li_values_random_dict[f'Li{n}'] = li_value

        else:
            self.li_values_dict[f'Li{n}'] = li_value
                
        return li_value
    
    def evaluate_theta_sharp_up_to_halfweight(self, halfweight, random_evaluation=False, eval_dict=None):
        '''
        Evaluate all values of log(), Li(n) up to n=halfweight. Evaluating these functions
        stores their outputs in a dictionary.
        '''
        
        self.log(random_evaluation, eval_dict)
        
        for n in range(1, halfweight + 1):
            if n > self.halfweight: break
            self.Li(n, random_evaluation, eval_dict)
            
    @staticmethod
    def compute_monomials_in_polylogs(halfweight, degree):
        '''
        Compute all monomials in the variables 'log', 'Li1', ...., 'Li{halfweight}' of weighted
        degree 'weight'.
        
        Return the output as a list of monomials, where each monomial is represented as a list of tuples
        of the form (variable, degree)
        '''
        variables = ['log'] + [f'Li{n}' for n in range(1, halfweight + 1)]
        degrees = [1] + [n for n in range(1, halfweight + 1) ]
        
        weighted_vectors = list(WeightedIntegerVectors(degree, degrees))
        
        monomials = []
        
        for degrees in weighted_vectors:
            monomial = []
            
            for index, variable in enumerate(variables):
                if degrees[index] == 0: continue
                
                monomial.append((variable, degrees[index]))
            
            monomials.append(monomial)
            
        return monomials

    def theta_sharp_matrix_in_given_degree(self, halfweight, degree, random_evaluation=False, eval_dict=None):
        '''
        Compute the matrix representing the action of theta^sharp the graded piece of O(\Pi_halfweight) of degree
        'degree', with respect to a monomial basis on O(\Pi_halfweight) and on O(U_S)[\Phi].
        
        The output will be a matrix whose coefficients are elements of O(U_S)_{\le halfweight}.
        
        If random_evaluation is true, the map is from O(\Pi_halfweight) to Q(\Phi), and the output will be
        an integer matrix (due to previous rescaling to ensure all coefficients are integral)
        
        '''
        self.evaluate_theta_sharp_up_to_halfweight(halfweight, random_evaluation, eval_dict)
        
        polylog_monomials = ThetaSharpOperator.compute_monomials_in_polylogs(halfweight, degree)
        phi_monomials = self.phi_algebra.monomials_of_degree(degree)
        phi_monomials.reverse()
        
        theta_sharp_values = []
        
        if random_evaluation:
            li_values_dict = self.li_values_random_dict
        
        else: 
            li_values_dict = self.li_values_dict
        
        #Evaluate the theta^sharp operator on each monomial in polylog_monomials
        for monomial in polylog_monomials:
            theta_sharp_val = 1
            
            for variable, degree in monomial:
                theta_sharp_val *= li_values_dict[variable]^degree
        
            theta_sharp_values.append(theta_sharp_val) 
            
        matrix = [[theta_sharp_val[phi_monomial] for phi_monomial in phi_monomials]
                  for theta_sharp_val in theta_sharp_values]
            
        return matrix
    
    def upper_bound_on_dimension_of_kernel(self, halfweight, degree, test_integers=False, clear=True):
        '''
        Compute a upper bound for the dimension of the kernel of theta^sharp, when restricted
        to polylogarithms up to Li_d and elements of \O(\Pi_d) of degree at most degree.
        
        While the output of this function is provably a lower bound, it is almost certainly
        an upper bound two, with the likelihood of success increasing on successive iterations
        '''
        
        if test_integers:
            self.random_eval_dict = self.set_random_integers(clear=True)
            
        else:
            self.random_eval_dict = self.generate_random_integers(clear=True)
        
        # If the halfweight is odd, then theta^sharp(Li_halfweight) contains \Phi_{\sigma_{halfweight}}, and that
        # is the only place that \Phi variable occurs. Hence, Li_halfweight cannot occur in any
        # Coleman function. So we can compute the matrix in one halfweight lower.
        
        if halfweight % 2 == 1:
            halfweight -= 1
        
        matrix = self.theta_sharp_matrix_in_given_degree(halfweight, degree, True, self.random_eval_dict)
        
        return Matrix(ZZ, matrix).left_nullity()

def main_parallel(depth, nprocesses):
    manager = Manager()
    cache = manager.dict()

    with Pool(processes=nprocesses) as pool:
        S = ThetaSharpOperator(depth, 2)
        shuf = S.shuffle_algebra
        words_to_compute = S.generate_words_to_precompute(depth)
 
        results = [pool.apply_async(shuf.to_dual_pbw_element_async, args=(word,cache)) for word in words_to_compute]
        results_dict = {r.get()[0] : r.get()[1] for r in results}

    shuf.cache = results_dict
    shuf.precomputed = True

    result = S.upper_bound_on_dimension_of_kernel(depth, 17, True, True)
    result2 = S.upper_bound_on_dimension_of_kernel(depth, 18, True, False)
    print(f"The number of Coleman functions in depth {depth} and weight 17 is {result}")
    print(f"The number of Coleman functions in depth {depth} and weight 18 is {result2}")

def main_non_parallel(depth):
    S = ThetaSharpOperator(depth, 2)
    result = S.upper_bound_on_dimension_of_kernel(depth, 17, True, True)
    result2 = S.upper_bound_on_dimension_of_kernel(depth, 18, True, False)
    print(f"The number of Coleman functions in depth {depth} and weight 17 is {result}")
    print(f"The number of Coleman functions in depth {depth} and weight 18 is {result2}")

def compute_upper_bound_on_dimension_of_kernel(depth, parallel=False, nprocesses=8):
    if paralell:
        main_parallel(depth, nprocesses)
    else:
        main_non_parallel(depth)
    
if __name__ == "__main__":
    compute_upper_bound_on_dimension_of_kernel(6)
    #compute_upper_bound_on_dimension_of_kernel(16, True, 95)
