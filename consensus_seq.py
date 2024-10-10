import re
import time


class DNA:
    def __init__(self, seq):
        self.seq = seq
        self.seq_len = len(seq)

    def all_kmer(self, k):
        self.kmers = []
        for i in range(self.seq_len):
            if i + k <= self.seq_len:
                self.kmers.append(self.seq[i : k + i])
        return self.kmers

    def unique_kmers(self, kmers=None):
        """Find unique kmers, to find kmers run all_kmer(k) first!"""
        if kmers:
            return set(kmers)
        else:
            return set(self.kmers)


class ShortestCommonSuperstring(DNA):
    A = []

    def __init__(self, kmers=None):
        if kmers != None:
            self.kmers = kmers
            self.A = [""] * len(kmers)
            self.scs()

        else:
            print(
                "Warning! No kmers provided. You can load sequences using load_seq() function."
            )

    def load_seq(self, seq, k):
        self.k = k
        self.seq = seq
        super().__init__(seq)
        self.kmers = super().all_kmer(k)
        # self.scs()

    def overlap(self, str1, str2):
        """Finding overlap in suffix of str1 and prefix of str1"""
        len1 = len(str1)
        len2 = len(str2)
        combined = ""
        for i in range(len1):
            if re.match(str1[i:], str2):
                # print(str1, str2, re.match(str1[i:], str2).span())
                offset = re.match(str1[i:], str2).span()[1]
                combined = str1 + str2[offset:]
                return [offset, combined]
                break

    def graph(self):
        self.uniq_kmers = super().unique_kmers(kmers=self.kmers)  # vertices
        # print("uniq_kmers: ", self.uniq_kmers)
        self.edges = []  #
        tmp = []
        for kmer1 in self.uniq_kmers:
            for kmer2 in self.uniq_kmers:
                if kmer1 != kmer2:
                    tmp = self.overlap(kmer1, kmer2)
                    if not tmp:
                        continue
                    elif tmp[0] == 0:
                        continue
                    else:
                        self.edges.append(
                            (kmer1, kmer2, tmp[0], tmp[1])
                        )  ##combined=str1+str2[offset:]

        # print(self.uniq_kmers, self.edges)

    def scs(self):
        self.graph()
        self.max = 0
        # print("edges: ", self.edges)
        for edge in self.edges:
            # print(edge)
            if self.max < edge[2]:
                self.max = edge[2]
        i = len(self.kmers)
        tmp = []
        for index, edge in enumerate(self.edges):
            if edge[2] == self.max:
                tmp = edge
                # del self.edges[i]
                # print(self.kmers)
                self.kmers.remove(edge[0])
                self.kmers.remove(edge[1])
                self.kmers.append(tmp[3])
                # print(self.kmers)
                self.scs()
                break
        if len(self.kmers) >= i:
            print("____________________________________________________")
            print(self.kmers[0])
            print("Supersting comun mas corto: ", self.kmers)
            print("Longitud del superstring: ", len(self.kmers[0]))
            return self.kmers
        i = len(self.kmers)


sequences = [
    "ATCCGTTGAAGCCGCGGGC",
    "TTAACTCGAGG",
    "TTAAGTACTGCCCG",
    "ATCTGTGTCGGG",
    "CGACTCCCGACACA",
    "CACAGATCCGTTGAAGCCGCGGG",
    "CTCGAGTTAAGTA",
    "CGCGGGCAGTACTT",
]
sequences = [
    "CCTCGAGTTAA",
    "GCCCGCGGCTTCAACGGAT",
    "TTAAGTACTGCCCG",
    "ATCTGTGTCGGG",
    "AAGTACTGCCCGCG",
    "TGTGTCGGGAGTCG",
    "CTCGAGTTAAGTA",
    "CCCGCGGCTTCAACGGATCTGTG",
]
# SCS = ShortestCommonSuperstring(["TACC", "ACTAC", "CGGACT", "ACGGA"])
print("""
CCTCGAGTTAA
-CTCGAGTTAAGTA
-------TTAAGTACTGCCCG-
---------AAGTACTGCCCGCG
----------------GCCCGCGGCTTCAACGGAT
-----------------CCCGCGGCTTCAACGGATCTGTG
---------------------------------ATCTGTGTCGGG
------------------------------------TGTGTCGGGAGTCG""")
start_time = time.time()
SCS = ShortestCommonSuperstring(sequences)
end_time = time.time()
exec_time = end_time - start_time
print("Tiempo de ejecución: ", exec_time)
