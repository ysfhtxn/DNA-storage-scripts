import reedsolo as rs
import numpy as np
import math

class DNAReedSolomonCode:
    def __init__(self, nucleotide_bases=None, base=None, gfnc=3, necc=10):
        if nucleotide_bases is None:
            nucleotide_bases = ['A', 'C', 'G', 'T', 'M', 'K', 'R', 'Y']  # Default nucleotide bases
        self.nucleotide_bases = nucleotide_bases
        self.base = base if base else len(nucleotide_bases)  # Automatically derive base if not provided
        self.gfnc = gfnc  # Galois Field Number of Characters
        self.necc = necc  # Number of Error Correcting Codes
        self.init_reedsolo()

    def init_reedsolo(self):
        # Calculate a prime polynomial suitable for the number of bases
        prim = rs.find_prime_polys(c_exp=math.ceil(np.log2(self.base)) * self.gfnc, fast_primes=True, single=True)
        rs.init_tables(c_exp=math.ceil(np.log2(self.base)) * self.gfnc, prim=prim)

    def decimal_to_multi_base(self, n):
        res = 0
        n = str(n)
        for i in range(len(n)):
            res += int(n[i]) * (self.base ** (len(n) - i - 1))
        return res

    def multi_base_to_decimal(self, n):
        res = []
        while n > 0:
            res.append(n % self.base)
            n = n // self.base
        res.reverse()
        return ''.join(str(i) for i in res).zfill(math.ceil(np.log2(self.base)))

    def encode(self, dna):
        msg_to_encode = []
        for i in range(0, len(dna), math.ceil(np.log2(self.base))):
            multi_base_block = ''.join([str(self.nucleotide_bases.index(dna[i:i + math.ceil(np.log2(self.base))][j])) for j in range(math.ceil(np.log2(self.base)))])
            msg_to_encode.append(self.decimal_to_multi_base(multi_base_block))
        gen = rs.rs_generator_poly_all(len(msg_to_encode))
        mesecc = rs.rs_encode_msg(msg_to_encode, self.necc, gen=gen[self.necc])
        mesecc_str = ''.join([self.multi_base_to_decimal(i) for i in mesecc])
        rsed_dna = ''.join([self.nucleotide_bases[int(i)] for i in mesecc_str])
        return rsed_dna

    def decode(self, rsed_dna):
        msg_to_decode = []
        for i in range(0, len(rsed_dna), math.ceil(np.log2(self.base))):
            multi_base_block = ''.join([str(self.nucleotide_bases.index(rsed_dna[i:i + math.ceil(np.log2(self.base))][j])) for j in range(math.ceil(np.log2(self.base)))])
            msg_to_decode.append(self.decimal_to_multi_base(multi_base_block))
        rmes, recc, errata_pos = rs.rs_correct_msg(msg_to_decode, self.necc)
        decoded_str = ''.join([self.multi_base_to_decimal(i) for i in rmes])
        decoded_dna = ''.join([self.nucleotide_bases[int(i)] for i in decoded_str])
        return decoded_dna
    
    def introduce_mutation_rate(self, dna, mutation_rate):
        dna = list(dna)
        for i in range(len(dna)):
            if np.random.rand() < mutation_rate:
                dna[i] = np.random.choice(self.nucleotide_bases)
        return ''.join(dna)

def main():
    dnars = DNAReedSolomonCode(necc=12)
    
    dna_msg     = 'TTTCTGTTGGTGCTGATATTGCTAACAGGACCAGGCGAAGAACAGGACCAGGCGAAGCRKGYYKYARGMCRYMRGGKMTYKMGCCAKYCGGTGKMRTGGCCTRAYGCGRTTMAAMGKAYGCTCGTAACGGCTTKGAAATTYAGG'
    encoded_dna = dnars.encode(dna_msg)
    errored_dna = dnars.introduce_mutation_rate(encoded_dna, mutation_rate=0.01)
    decoded_dna = dnars.decode(errored_dna)
    
    print("Original DNA:", dna_msg)
    print("Encoded  DNA:", encoded_dna)
    print("Errored  DNA:", errored_dna)
    print("Decoded  DNA:", decoded_dna)
    
if __name__ == '__main__':
    main()