/***************************************************************************
 *  Description:
 *      Experimental read mapper to explore read mapping and alignment
 *      algorithms.
 *
 *  History: 
 *  Date        Name        Modification
 *  2023-03-26  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <xtend/file.h>
#include <xtend/mem.h>
#include <biolibc/fasta.h>
#include <biolibc/fastq.h>

int     align(const char *ref_file, const char *reads_file);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char    *ref_file, *reads_file;
    
    switch(argc)
    {
	case 3:
	    ref_file = argv[1];
	    reads_file = argv[2];
	    break;
	
	default:
	    usage(argv);
    }
    
    return align(ref_file, reads_file);
}


/***************************************************************************
 *  Description:
 *      Basic algorithm to align reads to a sequence.
 *      Just beginning: need to explore other
 *      approaches besides brute force to make mapping feasible for
 *      large genomes and transcriptomes.
 *
 *  History: 
 *  Date        Name        Modification
 *  2023-03-26  Jason Bacon Begin
 ***************************************************************************/

inline ssize_t  match(const char seq[], size_t seq_len,
		      const char read[], size_t read_len)

{
    // Allow whole read to be compared at end of sequence
    const char  *last_start = seq + seq_len - read_len;
    const char  *sp, *rp, *seq_start;
    
    for (seq_start = seq; seq_start < last_start; ++seq_start)
    {
	rp = read, sp = seq_start;
	// Need only check one for null byte if they're the same
	while ( (*rp == *sp) && (*rp != '\0') )
	    ++rp, ++sp;
	
	if ( (*rp == 0) || (*sp == 0) )
	    return seq_start - seq;
    }
    
    return -1;
}


int     align(const char *ref_file, const char *reads_file)

{
    FILE        *ref_fp, *reads_fp, *align_fp;
    bl_fasta_t  **sequences;
    bl_fastq_t  read;
    size_t      seq_count, seq, seq_len, read_count, read_len;
    ssize_t     offset;
    char        *read_ptr, *seq_ptr;
    
    if ( (ref_fp = xt_fopen(ref_file, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", ref_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( (reads_fp = xt_fopen(reads_file, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", reads_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( (align_fp = xt_fopen("alignments.txt", "w")) == NULL )
    {
	fprintf(stderr, "Cannot open aligments.txt: %s\n", strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( (sequences = xt_malloc(10000, sizeof(bl_fasta_t *))) == NULL )
    {
	fputs("Could not allocate sequences pointer array.\n", stderr);
	return EX_UNAVAILABLE;
    }
    seq_count = 0;
    
    // FIXME: Check success
    sequences[0] = xt_malloc(1, sizeof(bl_fasta_t));
    bl_fasta_init(sequences[0]);
    while ( bl_fasta_read(sequences[seq_count], ref_fp) == BL_READ_OK )
    {
	// puts(BL_FASTA_DESC(&temp_seq));
	++seq_count;
	if ( (sequences[seq_count] = xt_malloc(1, sizeof(bl_fasta_t))) == NULL )
	{
	    fprintf(stderr, 
		    "Could not allocate FASTA sequence #%zu.\n", seq_count);
	    return EX_UNAVAILABLE;
	}
    }
    fclose(ref_fp);
    printf("%zu sequences loaded.\n", seq_count);
    
    read_count = 0;
    bl_fastq_init(&read);
    
    // FIXME: Stopping at 200 reads for quick test and timing
    while ( bl_fastq_read(&read, reads_fp) == BL_READ_OK && read_count < 200 )
    {
	read_ptr = BL_FASTQ_SEQ(&read);
	for (seq = 0; seq < seq_count; ++seq)
	{
	    seq_ptr = BL_FASTA_SEQ(sequences[seq]);
	    seq_len = BL_FASTA_SEQ_LEN(sequences[seq]);
	    read_len = BL_FASTQ_SEQ_LEN(&read);
	    if ( (offset = match(seq_ptr, seq_len, read_ptr, read_len)) > 0 )
		fprintf(align_fp, "\ns[%zu], %zd\n", seq, offset);
	}
	if ( read_count % 100 == 0 )
	{
	    printf("%zu\r", read_count);
	    fflush(stdout);
	}
	++read_count;
    }
    fclose(reads_fp);
    fclose(align_fp);
    printf("%zu reads processed.\n", read_count);
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s reference.fa[.gz|.bz2|.xz] reads.fq[.gz|.bz2|.xz]\n", argv[0]);
    exit(EX_USAGE);
}
