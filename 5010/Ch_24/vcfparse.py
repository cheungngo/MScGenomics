"""Personal module for parsing VCF"""

import re
import io

class get_vcf:
    
    def __init__(self, filename):
        chr = dict()
        no = sum(1 for line in io.open(filename))
        with io.open(filename) as fh:
            for i in range(no):
                ln = fh.readline().strip()
                if not (re.search(r"^#", ln)):
                    ln = ln.split("\t")
                    if ln[0] not in chr.keys():
                        chr[ln[0]] = chromosome(ln[0])
                    chr[ln[0]].add_snp(ln[0],ln[1],ln[2],ln[3],ln[4])
        self.chr = chr

class SNP:
    
    def __init__(self, chrom, pos, id, ref, alt):
        assert ref != alt, "Error: ref == alt"
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        
    def is_transition(self):
        if self.ref == "A" or self.ref == "G":
            if self.alt == "A" or self.alt == "G":
                return(True)
        if self.ref == "C" or self.ref == "T":
            if self.alt == "C" or self.alt == "T":
                return(True)
        return(False)
    
    def is_transversion(self):
        if self.is_transition():
            return(False)
        return(True)
    
class chromosome:
    
    def __init__(self, chrom):
        self.chrom = chrom
        self.snploc = dict()
        
    def add_snp(self, chrom, pos, id, ref, alt):
        assert pos not in self.snploc.keys(), "Error: duplicated SNP"
        assert chrom == self.chrom, "Error: wrong chromosome"
        newsnp = SNP(chrom, pos, id, ref, alt)
        self.snploc[pos] = newsnp
        
    def count_transition(self):
        count = 0
        locs = self.snploc.keys()
        for loc in locs:
            if self.snploc[loc].is_transition():
                count += 1
        return(count)
    
    def count_transversion(self):
        count = 0
        locs = self.snploc.keys()
        for loc in locs:
            if self.snploc[loc].is_transversion():
                count += 1
        return(count)
    
    def density(self, l, m):
        count = 0
        for loc in self.snploc.keys():
            loc = int(loc)
            if loc >= l and loc <= m:
                count += 1
        den = count / (m - l + 1) * 1000
        return(den)
    
    def best_region(self, region_size):
        
        lastsnp = int(sorted(self.snploc.keys(), reverse = True)[0])
        best = [0.0, 1, region_size - 1] # density, start, end
        
        for loc in range(1, lastsnp, region_size):
            loc = int(loc)
            den = self.density(loc, loc + region_size - 1)
            if den > best[0]:
                best = [den, loc, loc + region_size - 1]
        return(best)

    def region_trans(self, l, m):
        count_snp = 0
        count_transition = 0
        count_transversion = 0
        for loc in self.snploc.keys():
            loc = int(loc)
            if loc >= l and loc <= m:
                count_snp += 1
                if self.snploc[str(loc)].is_transition():
                    count_transition += 1
                if self.snploc[str(loc)].is_transversion():
                    count_transversion += 1
        try:
            return [count_snp, round(count_transition / count_snp,2), round(count_transversion / count_snp,2)]
        except:
            return [count_snp, 0, 0]

    def region_den(self, region_size):
        
        lastsnp = int(sorted(self.snploc.keys(), reverse = True)[0])
        rden = dict()
        
        for loc in range(1, lastsnp, region_size):
            loc = int(loc)
            den = round(self.density(loc, loc + region_size - 1),2)
            reg_trans = self.region_trans(loc, loc + region_size - 1)
            no_snp = reg_trans[0]
            percent_transition = reg_trans[1]
            key = str(loc) + ".." + str(loc + region_size - 1)
            rden[key] = [key, den, percent_transition, no_snp]
        return(rden)