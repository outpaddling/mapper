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
 *      Basic algorithm to align reads to a sequence.  More than twice
 *      as fast as strcmp().  Just beginning: need to explore other
 *      approaches besides brute force to make mapping feasible for
 *      large genomes and transcriptomes.
 *
 *  History: 
 *  Date        Name        Modification
 *  2023-03-26  Jason Bacon Begin
 ***************************************************************************/

inline int  match(const char * restrict p1, const char * restrict p2)

{
    while ( (*p1 == *p2) && (*p1 != '\0') )
	++p1, ++p2;
    
    return *p2 - *p1;
}


int     align(const char *ref_file, const char *reads_file)

{
    FILE        *ref_fp, *reads_fp;
    bl_fasta_t  temp_seq, **sequences;
    bl_fastq_t  read;
    size_t      seq_count, seq, read_count;
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
    
    if ( (sequences = xt_malloc(10000, sizeof(bl_fasta_t *))) == NULL )
    {
	fputs("Could not allocate sequences pointer array.\n", stderr);
	return EX_UNAVAILABLE;
    }
    seq_count = 0;
    
    bl_fasta_init(&temp_seq);
    while ( bl_fasta_read(&temp_seq, ref_fp) == BL_READ_OK )
    {
	// puts(BL_FASTA_DESC(&temp_seq));
	if ( (sequences[seq_count] = xt_malloc(1, sizeof(bl_fasta_t))) == NULL )
	{
	    fprintf(stderr, 
		    "Could not allocate FASTA sequence #%zu.\n", seq_count);
	    return EX_UNAVAILABLE;
	}
	*sequences[seq_count] = temp_seq;
	++seq_count;
    }
    fclose(ref_fp);
    printf("%zu sequences loaded.\n", seq_count);
    
    read_count = 0;
    bl_fastq_init(&read);
    while ( bl_fastq_read(&read, reads_fp) == BL_READ_OK )
    {
	read_ptr = BL_FASTQ_SEQ(&read);
	for (seq = 0; seq < seq_count; ++seq)
	{
	    for (seq_ptr = BL_FASTA_SEQ(sequences[seq]); *seq_ptr != 0; ++seq_ptr)
	    {
		// Checking first char before calling match() speeds it up a bit
		if ( *read_ptr == *seq_ptr )
		    if ( match(read_ptr, seq_ptr) == 0 )
			printf("\nMatch found in sequence %zu.\n", seq);
	    }
	}
	if ( read_count % 100 == 0 )
	{
	    printf("%zu\r", read_count);
	    fflush(stdout);
	}
	++read_count;
    }
    fclose(reads_fp);
    printf("%zu reads processed.\n", read_count);
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n", argv[0]);
    exit(EX_USAGE);
}
