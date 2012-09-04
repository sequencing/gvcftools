
gvcftools - utilities for gVCF files


Chris Saunders (csaunders@illumina.com)
Version: ${VERSION}


SUMMARY:

gvcftools is a set of utilities applicable to VCF files following the
gVCF convention for representation of non-variant positions


Currently available tools are:

trio - used to count and report inheritance conflicts and joint
       coverage in three gVCF files representing a parent-child trio

twins - used to count genotype conflicts and joint coverage between
        two samples, typically used for technical replicates or
        monozygotic twins


INSTALLATION:

The tarball includes the basic program dependendencies, so all that
should be required is is to extract the tarball, move to the root
package directory and run 'make'. Final binaries will be installed in
'$(PACKAGE_DIR)/bin'

Note that the dependendencies included are subsets of the tabix,
samtools and boost packages.


