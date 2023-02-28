from __future__ import print_function
import numpy as np
from fractions import Fraction
import itertools
from collections import Counter

def expand(term):
#Recursive function that takes as input a term and outputs a list of terms
#By term we mean a list consisting of a Fraction followed by a nonzero number of lists of strings
    for i in range(len(term)-1,0,-1):   
    #Inspect the elements of term from right to left, excluding term[0]
        if term[i][0]=='bD':
        #If term[i] is of the form ['bD','j'], where j is any index
            if i==len(term)-1:
            #If this 'bD' is the last element of term
                return expand(term[0:i]+[[','],['D',term[i][1]],['x']])
            else:
            #If not, then we assume it is followed by the element [','] 
            #and then at least one more element (so i+2<len(term)).
            #We produce three terms which are in turn expanded.
                if [','] in term[i+2:]:
                #If there is another comma after the comma following this 'bD'
                    index_next_comma = term[i+2:].index([','])+i+2
                else:
                    index_next_comma = len(term)

                return(
                    expand(term[0:i]+[[','],['D',term[i][1]],['x']]+term[i+1:])
                    +expand(term[0:i]+[[','],['D',term[i][1]]]+term[i+2:])
                    +expand(term[0:i]+term[i+1:index_next_comma]+[term[i]]+term[index_next_comma:])
                    )                
    #If there is no 'bD' in term, we end up here and we just return [term]
    return [term]

def term_to_string(term,k):
    #example:   term = [Fraction(1,3),['a'],['D','i'],['x'],[','],['b']]
    #           k = 4
    #returns:   '(1/3)T^x_{F_{4,d}^{[2]}}(aD_ix,b)'
    string = ''
    for f in term:
        if type(f) is Fraction:
            if f is not Fraction(1,1) and f!=1:
                string+='('+str(f)+')'
            string+='T^x_{F_{'+str(k)+',d}^{['+str(term.count([','])+1)+']}}('
        elif len(f)>1:
            string+=str(f[0])+'_'+str(f[1])
        else:
            string+=str(f[0])
    string+=')'
    return string

def sort_indices(string,alphabet):
    #example:   string = 'T^x_{f^{[4]}}(xD_jD_kx,a_i,D_ja_l,a_lD_kD_jx)'
    #           alphabet = 'ijkl'
    #returns:   'T^x_{f^{[4]}}(xD_iD_jx,a_k,D_ia_l,a_lD_jD_ix)'
    count = 0
    for i in np.arange(len(string)):
        if string[i] in alphabet[count:]:
            if string[i] == alphabet[count]:
                count += 1
            else:
                char = string[i]
                string = string.replace(alphabet[count],'@') #Swap all occurences of alphabet[count] and char.
                string = string.replace(char,alphabet[count])
                string = string.replace('@',char) # @ is a placeholder
                count += 1
            if count == len(alphabet):
                return string
    return string

def V(index):
    return [    [['x'],['bD',index]],
                [['x'],['bD',index]],
                [['a',index]]           ]

def P(index):
    return [    [['x'],['bD',index],['bD',index]],
                [['a',index],['bD',index]],
                [['a/d']]                           ]

def print_I(k):
#Function which computes and prints an explicit expression for I_k of the NC d-torus,
#leaving d as an indetermined variable. k is assumed a nonnegative integer
    if k%2==1:
        print('$$I_'+str(k)+'=0$$')
        return 0
    if k==0:
        print('$$\pi^{-d/2}I_0=x^{-d/2}$$')
        return 1
    if k%4==2:
        print('$$-\pi^{-d/2}I_'+str(k)+'=$$')
    else:
        print('$$\pi^{-d/2}I_'+str(k)+'=$$')
    k_over_2 = int(k/2)
    how_many_terms = 0
    for m in range(k_over_2,k+1):
        if m>k_over_2:
            print('$$+',end='')
        else:
            print('$$',end='')
        indices_all = [chr(105+i) for i in range(m)]
        numel_A = 2*m-k
        if numel_A == 2:
            indices_inA = [chr(105),chr(105)]
            indices_all.remove(chr(106))
            print('\sum_{',','.join(indices_all),'}\\frac{1}{d}\Bigg(',sep='')
        else:
            indices_inA = [chr(105+i) for i in range(numel_A)]
            if numel_A == 0:
                print('\sum_{',','.join(indices_all),'}\Bigg(',sep='')
            else:
                print('\sum_{',','.join(indices_all),'}c_d^{(',','.join(indices_inA),')}\Bigg(',sep='')
        indices_notinA = [chr(105+i) for i in range(numel_A,m)]
        terms_m = []
        for A in itertools.combinations(range(1,m+1),2*m-k):
            w_A = []
            ind_counter_inA = 0
            ind_counter_notinA = 0
            for i in range(m):
                if i+1 in A:
                    w_A.append(V(indices_inA[ind_counter_inA]))
                    ind_counter_inA += 1
                else:
                    w_A.append(P(indices_notinA[ind_counter_notinA]))
                    ind_counter_notinA += 1
            for p in itertools.product([0,1,2],repeat=m):
                term = [Fraction(1,1)]
                for i in range(m):
                    if i>0:
                        term.append([','])
                    term.extend(w_A[i][p[i]])
                terms_m.extend(expand(term))
                
        #Make strings, sort indices, put in the list termstrings_m
        termstrings_m = []
        for term in terms_m:
            termstrings_m.append(sort_indices(sort_indices(term_to_string(term,k),indices_inA),indices_notinA))
                
        keys = list(Counter(termstrings_m).keys())
        values = list(Counter(termstrings_m).values())
        how_many_terms += len(values)
        for i in range(len(values)):
            if i>0:
                print('+',end='')
            if values[i] != 1:
                print(values[i],end='')
            print(keys[i],end='')
            if i%3==2 and i!=(len(values)-1):
                print('$$')
                print('$$',end='')
        print('\Bigg)$$')
    return how_many_terms

#Some examples of application of the above code:          
#example_term = [Fraction(1,3),['a'],['bD','j'],['bD','l'],[','],['D','i'],['x'],[','],['b']]
#example_term_2 = [Fraction(1,1),['a'],['bD','i'],['bD','j'],[','],['b'],['bD','k'],[','],['c']]
#print(len(expand(example_term_2)))
#how_many_terms_in_I_k = print_I(2)
#print('Amount of terms = '+str(how_many_terms_in_I_k))

print_I(4)
