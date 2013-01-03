#!/usr/bin/env python
from pykent.common.Sanity import errAbort, canBeInt, canBeNum, roughlyEqual
from pykent.common.DNAUtil import isDNA

MuTectorVersion = "## muTector v1.0.27200"
MuTectorColumns = "contig position ref_allele alt_allele tumor_name normal_name score dbsnp_site covered power tumor_power normal_power total_pairs improper_pairs map_Q0_reads t_lod_fstar tumor_f contaminant_fraction contaminant_lod t_ref_count t_alt_count t_ref_sum t_alt_sum t_ins_count t_del_count normal_best_gt init_n_lod n_ref_count n_alt_count n_ref_sum n_alt_sum judgement"
NumMuTectorCols = len(MuTectorColumns.split(' '))
MuTectorAlgorithm = "MuTect"

class MuTectorLine:

	def __init__(self, line, snpAnnotationFilters):
		''' Initialize a MuTectorLine object with all columns properly instantiated.
			Assumes columns are in the exact ordering given by the Columns variable. '''
		if line.endswith('\n'):
			line = line[:-1]
		line = line.split('\t')
		if len(line) != NumMuTectorCols:
			errAbort("MuTectorLine has %d columns (not %d): %s" % (len(line), self.NumCols, "\t".join(line)))

		self.contig = line[0]          # Chromosome
		self.position = int(line[1])   # Chrom position
		self.ref_allele = line[2]      # DNA
		self.alt_allele = line[3]      # DNA
		self.tumor_name = line[4]
		self.normal_name = line[5]
		self.score = int(line[6])
		self.dbsnp_site = line[7]
		self.covered = line[8]
		self.power = float(line[9])
		self.tumor_power = float(line[10])
		self.normal_power = float(line[11])
		self.total_pairs = int(line[12])
		self.improper_pairs = int(line[13])
		self.map_Q0_reads = int(line[14])
		self.t_lod_fstar = float(line[15])
		self.tumor_f = float(line[16])
		self.contaminant_fraction = float(line[17])
		self.contaminant_lod = float(line[18])
		self.t_ref_count = int(line[19])
		self.t_alt_count = int(line[20])
		self.t_ref_sum = int(line[21])
		self.t_alt_sum = int(line[22])
		self.t_ins_count = int(line[23])
		self.t_del_count = int(line[24])
		self.normal_best_gt = line[25]   # DNA
		self.init_n_lod = float(line[26])
		self.n_ref_count = int(line[27])
		self.n_alt_count = int(line[28])
		self.n_ref_sum = int(line[29])
		self.n_alt_sum = int(line[30])
		self.judgement = line[31]
		self.algo = MuTectorAlgorithm
		self.validateMuTectorLine()
		self.makeOurJudgment(snpAnnotationFilters)

	def validateMuTectorLine(self):
		''' Make sure all MuTector fields make sense '''
		if self.position < 0:
			errAbort("Position must be non-negative: %s" % self.__str__())
		if not isDNA(self.ref_allele):
			errAbort("Ref allele must be DNA: %s" % self.__str__()) 
		if not isDNA(self.alt_allele):
			errAbort("Alt allele must be DNA: %s" % self.__str__()) 
		if self.score < 0:
			errAbort("Score must be non-negative: %s" % self.__str__())
		if self.power < 0 or self.power > 1:
			errAbort("Power must be in [0,1]: %s" % self.__str__())
		if self.tumor_power < 0 or self.tumor_power > 1:
			errAbort("Tumor power must be in [0,1]: %s" % self.__str__())
		if self.normal_power < 0 or self.normal_power > 1:
			errAbort("Normal power must be in [0,1]: %s" % self.__str__())
		if self.total_pairs < 0:
			errAbort("Total pairs must be non-negative: %s" % self.__str__())
		if self.improper_pairs < 0 or self.improper_pairs > self.total_pairs:
			errAbort("Improper pairs must be non-negative and less than total pairs: %s" % self.__str__())
		if self.map_Q0_reads < 0:
			errAbort("map_Q0_reads must be non-negative: %s" % self.__str__())
		if self.contaminant_fraction < 0 or self.contaminant_fraction > 1:
			errAbort("Contaminant fraction must be in [0,1]: %s" % self.__str__())
		if self.t_ref_count < 0 or self.t_ref_sum < 0 or self.t_ref_count > self.t_ref_sum:
			errAbort("Tumor ref count/sum are bad: %s" % self.__str__())
		if self.t_alt_count < 0 or self.t_alt_sum < 0 or self.t_alt_count > self.t_alt_sum:
			errAbort("Tumor alt count/sum are bad: %s" % self.__str__())
		if self.t_ins_count < 0 or self.t_del_count < 0:
			errAbort("Tumor insertion/deletion counts must be non-negative: %s" % self.__str__())
		if not isDNA(self.normal_best_gt):
			errAbort("Normal best GT must be DNA: %s" % self.__str__())
		if self.n_ref_count < 0 or self.n_ref_sum < 0 or self.n_ref_count > self.n_ref_sum:
			errAbort("Normal ref count/sum are bad: %s" % self.__str__())
		if self.n_alt_count < 0 or self.n_alt_sum < 0 or self.n_alt_count > self.n_alt_sum:
			errAbort("Normal alt count/sum are bad: %s" % self.__str__())

	def t_var_freq(self):
		''' Return the tumor variant frequency '''
		if self.t_alt_count+self.t_ref_count == 0:
			return 0
		return self.t_alt_count * 1./(self.t_alt_count + self.t_ref_count)

	def n_var_freq(self):
		''' Return the normal variant frequency '''
		if self.n_alt_count + self.n_ref_count == 0:
			return 0
		return self.n_alt_count * 1./(self.n_alt_count + self.n_ref_count)

	def getKnownVariantStatus(self):
		''' Return a string representing the known variant status of the call '''
		if self.dbsnp_site in ('DBSNP', 'UNCOVERED'):
			return "MuTect_" + self.dbsnp_site
		return "NOVEL"

	def makeOurJudgment(self, snpFilters):
		''' Determine whether we should label the mutation as yes/no/maybe based on the info we have and our filtering criteria '''
		self.ourJudgment = "yes"
		judgmentReasons = []
		if (self.t_ref_count + self.t_alt_count) < snpFilters['minTumorReads']:
			self.ourJudgment = "no"
			judgmentReasons.append("TUMOR_READ_COUNT_TOO_LOW")
		if self.t_alt_count < snpFilters['minTumorAltReads']:
			self.ourJudgment = "no"
			judgmentReasons.append("TUMOR_ALT_READ_COUNT_TOO_LOW")
		if self.t_var_freq() < snpFilters['minTumorVarFreq']:
			self.ourJudgment = "no"
			judgmentReasons.append("TUMOR_VAR_FREQ_TOO_LOW")
		if (self.n_ref_count + self.n_alt_count) < snpFilters['minNormalReads']:
			self.ourJudgment = "no"
			judgmentReasons.append("NORMAL_READ_COUNT_TOO_LOW")
		if self.n_alt_count > snpFilters['maxNormalAltReads']:
			self.ourJudgment = "no"
			judgmentReasons.append("NORMAL_ALT_READ_COUNT_TOO_HIGH")
		if self.n_var_freq() > snpFilters['maxNormalVarFreq']:
			self.ourJudgment = "no"
			judgmentReasons.append("NORMAL_VAR_FREQ_TOO_HIGH")
		if self.n_var_freq() >= self.t_var_freq():
			self.ourJudgment = "no"
			judgmentReasons.append("NORMAL_VAR_FREQ_HIGHER_THAN_TUMOR_VAR_FREQ")
		if self.covered not in snpFilters['coveredAllowed']:
			self.ourJudgment = "no"
                        judgmentReasons.append("POSITION_NOT_COVERED")

		self.ourJudgmentReasons = ",".join(judgmentReasons)

	def shouldBeKept(self, snpRemovalFilters):
		''' Return True iff we should keep this mutation in our output file. '''
		return (self.judgement in snpRemovalFilters['muTectJudgmentsAllowed'])

	def AnnoVarStr(self):
		''' Write out a description of this mutation that is suitable for AnnoVar parsing '''
		return "%s\t%d\t%d\t%s\t%s" % (self.contig, self.position, self.position, self.ref_allele, self.alt_allele)		
	def __str__(self):
		''' Return a string representation of a MuTectorLine object '''
		return 

def validateMutectorFile(version, cols):
	''' Abort if MuTect file is not a recognized version or does not have the proper headers '''
	global MuTectorVersion, MuTectorColumns, NumMuTectorColumns
	if version.strip() != MuTectorVersion:
		print "Warning: MuTector version (%s) not what we are expecting (%s)..." % (version, MuTectorVersion)
	if cols != MuTectorColumns:
		print "MuTectorc olumns not the expected columns."
		actualCols = cols.split(' ')
		expectedCols = MuTectorColumns.split(' ')
		actualLen = len(actualCols)
		print "#Col\tActual\tExpected"
		for i in xrange(max(actualLen, NumMutectorColumns)):
			a = ""
			if i < actualLen:
				a = actualCols[i]
			e = ""
			if i < expectedLen:
				e = expectedCols[i] 
			note = ""
			if a != e:
				note = "*"
			print "%d\t%s\t%s\t%s" % (i, a, e, note)
		errAbort("")

