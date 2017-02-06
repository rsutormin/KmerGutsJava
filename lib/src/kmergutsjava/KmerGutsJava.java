package kmergutsjava;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringTokenizer;

/*

kmer_guts.c can be compiled into either a 5-mer or an 8-mer version.  I have 
labeled the critical changes with 

     ### CHANGE THIS FOR DIFFERING Ks ###

Note: to run use

   build_data_directory KmerData 8

Then use

   kmer_guts -w -D KmerData < input.contigs

to build an memory map -- it can take 10-15 minutes.

Then use

   kmer_guts -D KmerData < input.contigs

to run.  Eventually, the memory-mapped hash table, which is large,
will become resident, and these should go fairly quickly.

   kmer_guts -D KmerData -a  < amino_acid_sequences.fasta > hits

can be used to call translated PEGs.

where KmersData is a directory that must contain

       final.kmers     [a file of [kmer,avg-offset-from-end,function-index,otu-index]]
       function.index  [a file of [index,function] pairs]
       otu.index       [a file of [index,otu] pairs]
----------------

Conceptually, the data associated with each signature Kmer is

        1. the protein kmer
        2. the average offset from the end of the protein
        3. a set of numeric values that include

                function index
                OTU index

The program takes as input a  file that should be thought of as a sequence
of fasta files in which a single line

          >FLUSH

occurs after each infut file.  That is,

          >contig1
      acgtacgt
      >FLUSH
      >contig2
      acgtacgt

is an input file containing two input fasta files (each containing a single, short
contig).  The '>FLUSH' causes all files to be flushed.

In effect, the program processes a sequence of requests.  The output for each request
is a piece of stadout terminated by a line

          //

The lines in the output for a request are of four forms.  The first type of line begins
the processing of a contig and looks like

         processing contig[length]

Then, there will be a message like

         TRANSLATION contig length-contig frame offset-of-start  

before beginning output for each of the six frames.

Thus, there will be six such lines printed for each input contig.  Following
each of these "TRANSLATION..." lines, there will be a set of lines like

         CALL Start-in-protein-seq end-in-protein-seq number-hits function-index function

Finally, after all six frames have been processed (for DNA - just one sequence for aa input),
you get a line of the form

         OTU-COUNTS cnt1-otu1 cnt2-otu2 ...

This code uses a table indicating which K-mers are signatures. I call this the 
"kmer_bits" table.  You can run this code for K ranging from 5 to 8.

#################
The line

    #define K 8

defines the kmer size.  
The code is intended to work with only 5-mers or 8-mers.  That is:

YOU MUST SET K TO EITHER 5 OR 8 (sorry).
############################################

COMMAND LINE ARGUMENTS:

    -a         means amino acid input sequence (defaults to DNA)

    -d Level   sets debugging level (1 shows hits; after that it gets intense

    -m MinHits minimum number of hits to get CALLed

    -M MinWeighted  Kmers are now weighted; this is the min of the sum of the weights,

    -O         use order contraint

    -g MaxGap  sets maximum allowed gap between HITS

    -D Data    sets the Data directory where the memory map lives

    -s HashSize make sure that the value is the same when you save the memory map
                and when you use it to search

    -w          write the memory map (means Data must contain final.kmers and the indexes

 */
public class KmerGutsJava {

//    #define _FILE_OFFSET_BITS 64
//    #include <stdio.h>
//    #include <stdlib.h>
//    #include <sys/socket.h>
//    #include <netinet/in.h>
//    #include <arpa/inet.h>
//    #include <ctype.h>
//    #include <fcntl.h>
//    #include <unistd.h>
//    #include <signal.h>
//    #include <sys/types.h>
//    #include <sys/mman.h>
//    #include <sys/stat.h>
//    #include <errno.h>
//    #include <string.h>

    /* parameters to main -- accessed globally */
    static int debug = 0;
    static int aa    = 0;
    static long size_hash =  1400303159; /* 1400303159  tot_lookups=13474100 retry=2981020 for 5.contigs 4.684 sec */
    /* 2147483648  tot_lookups=13474100 retry=1736650  */
    /* 1073741824  tot_lookups=13474100 retry=4728020  */
    static int write_mem_map = 0;
    static String data_dir;

    public static final int K = 8;
    public static final int MAX_SEQ_LEN = 500000000;

    public static final long CORE = 20L*20L*20L*20L*20L*20L*20L;

    public static final long MAX_ENCODED = CORE*20L; 

    static int tot_lookups = 0;
    static int retry  = 0;


    public static final char[] genetic_code = {
            'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
            'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
            'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
            '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
    };

    public static final char[] prot_alpha = {
            'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
    };


    public static class sig_kmer {
        long which_kmer;
        int  otu_index;
        int  avg_from_end;
        int  function_index;
        float function_wt;
    } // sig_kmer_t;

    public static final int VERSION = 1;
    
    public static class kmer_memory_image {
        long num_sigs;
        long entry_size;
        long  version;
    } // kmer_memory_image_t;

    public static class kmer_handle {
        sig_kmer[] kmer_table;
        long num_sigs;
        String[] function_array;   /* indexed by fI */
        String[] otu_array;        /* OTU indexes point at a representation of multiple OTUs */
    } // kmer_handle_t;

    /* the following stuff was added to condense sets of hits to specific calls.
       The basic strategy is to use a set of global variables to retain state and flush
       calls (bad, bad, bad...).
     */

    public static class hit {
        int   oI;
        int from0_in_prot;      /* offset from start of protein sequence */
        short avg_off_from_end;  /* average offset from the end */
        int fI;
        float function_wt;
    } // hit_t;

    public static final int MAX_HITS_PER_SEQ = 40000;
    static hit[] hits = new hit[MAX_HITS_PER_SEQ]; 
    static int   num_hits = 0;

    public static final int OI_BUFSZ = 5;
    
    public static class otu_count {
        int oI;
        int count;
    } 
    
    otu_count[] oI_counts = new otu_count[OI_BUFSZ];
    static int num_oI = 0;

    static int   current_fI;
    static char[]  current_id = new char[300];
    static int   current_length_contig;
    static char  current_strand;
    static int   current_prot_off;
    static int   order_constraint = 0;
    static int   min_hits = 5;
    static int   min_weighted_hits = 0;
    static int   max_gap  = 200;

    public static final int MAX_FUNC_OI_INDEX = 1000000;
    public static final int MAX_FUNC_OI_VALS  = 100000000;

    /* =========================== end of reduction global variables ================= */

    public static byte to_amino_acid_off(char c) {
        switch (c)
        {
        case 'A':
            return 0;

        case 'C':
            return 1;

        case 'D':
            return 2;

        case 'E':
            return 3;

        case 'F':
            return 4;

        case 'G':
            return 5;

        case 'H':
            return 6;

        case 'I':
            return 7;

        case 'K':
            return 8;

        case 'L':
            return 9;

        case 'M':
            return 10;

        case 'N':
            return 11;

        case 'P':
            return 12;

        case 'Q':
            return 13;

        case 'R':
            return 14;

        case 'S':
            return 15;

        case 'T':
            return 16;

        case 'V':
            return 17;

        case 'W':
            return 18;

        case 'Y':
            return 19;
        }
        return 20;
    }

    public static char compl(char c)
    {
        switch (c)
        {
        case 'a':
            return 't';
        case 'A':
            return 'T';

        case 'c':
            return 'g';
        case 'C':
            return 'G';

        case 'g':
            return 'c';
        case 'G':
            return 'C';

        case 't':
        case 'u':
            return 'a';
        case 'T':
        case 'U':
            return 'A';

        case 'm':
            return 'k';
        case 'M':
            return 'K';

        case 'r':
            return 'y';
        case 'R':
            return 'Y';

        case 'w':
            return 'w';
        case 'W':
            return 'W';

        case 's':
            return 'S';
        case 'S':
            return 'S';

        case 'y':
            return 'r';
        case 'Y':
            return 'R';

        case 'k':
            return 'm';
        case 'K':
            return 'M';

        case 'b':
            return 'v';
        case 'B':
            return 'V';

        case 'd':
            return 'h';
        case 'D':
            return 'H';

        case 'h':
            return 'd';
        case 'H':
            return 'D';

        case 'v':
            return 'b';
        case 'V':
            return 'B';

        case 'n':
            return 'n';
        case 'N':
            return 'N';

        }
        return c;
    }


    public static char[] rev_comp(char[] data) {
        int n = data.length;
        char[] cdata = new char[n];
        int p  = n - 1;
        int pc = 0;
        while (n-- > 0) {
            cdata[pc++] = compl(data[p--]);
        }
        return cdata;
    }

    public static long encoded_kmer(byte[] data, int pos) {
        long encodedK = 0;
        for (int i = 0; i < K; i++) {
            int add = data[pos + i];
            if (add >= 20) {
                return -1;
            }
            encodedK = (encodedK * 20) + add;
        }
        if (encodedK > MAX_ENCODED) {
            String error = "bad encoding (" + encodedK + " > " + MAX_ENCODED + ") - " +
            		"input must have included invalid characters:";
            for (int i=0; i < K; i++) {
                error += " " + data[pos + i];
            }
            throw new IllegalStateException(error);
        }
        return encodedK;
    }

    public static int dna_char(char c)
    {
        switch (c)
        {
        case 'a':
        case 'A':
            return 0;

        case 'c':
        case 'C':
            return 1;

        case 'g':
        case 'G':
            return 2;

        case 't':
        case 'u':
        case 'T':
        case 'U':
            return 3;

        }
        return 4;
    }

    public static void translate(char[] seq, int off, char[] pseq, byte[] pIseq) {

        int i;
        int max = seq.length - 3;
        int p = 0;
        for (i=off; (i <= max); ) {
            int c1 = dna_char(seq[i++]);
            int c2 = dna_char(seq[i++]);
            int c3 = dna_char(seq[i++]);
            if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
                int I = (c1 * 16) + (c2 * 4) + c3;
                char prot_c = genetic_code[I];
                pseq[p] = prot_c;
                pIseq[p] = to_amino_acid_off(prot_c);
            }
            else {
                pseq[p]  = 'x';
                pIseq[p] = 20;
            }
            p++;
        }
        if (p < pseq.length) {
            pseq[p] = 0;
            pIseq[p] = 21;
        }
    }

    /*List<String> load_indexed_ar(File filename, int[] retSize) {
        char **index_ar = malloc(MAX_FUNC_OI_INDEX * sizeof(char *));
        char *vals      = malloc(MAX_FUNC_OI_VALS);
        char *p         = vals;
        FILE *ifp      = fopen(filename,"r");
        if (ifp == NULL) { 
            fprintf(stderr,"could not open %s\n",filename);
            exit(1);
        }

        int sz = 0;
        int j;
        while ((fscanf(ifp,"%d\t",&j) == 1) && fgets(p,1000,ifp)) {
            if (sz != j) {
                fprintf(stderr,"Your index must be dense and in order (see line %ld, should be %d)\n",p-vals,sz);
                exit(1);
            }
            index_ar[sz] = p;
            p += strlen(index_ar[sz]) -1;
            *(p++) = '\0';
            if ((sz >= MAX_FUNC_OI_INDEX) || ((p-vals) > (MAX_FUNC_OI_VALS - 1000))) {
                fprintf(stderr,"Your function or oI index arrays are too small; bump MAX_FUNC_OI_INDEX and MAX_FUNC_OI_VALS\n");
                exit(1);
            }

            sz += 1;
        }
        if (retSize != null) {
            retSize[0] = sz;
        }
        return index_ar;
    }

    List<String> load_functions(File file) {
        return load_indexed_ar(file, null);
    }

    List<String> load_otus(File file) {
        return load_indexed_ar(file, null);
    }*/

    /*long lookup_hash_entry(sig_kmer[] sig_kmers, long encodedK) {
        long long  hash_entry = encodedK % size_hash;
        if (debug >= 2)
            tot_lookups++;
        while ((sig_kmers[hash_entry].which_kmer <= MAX_ENCODED) && (sig_kmers[hash_entry].which_kmer != encodedK)) {
            if (debug >= 2)
                retry++;
            hash_entry++;
            if (hash_entry == size_hash)
                hash_entry = 0;
        }
        if (sig_kmers[hash_entry].which_kmer > MAX_ENCODED) {
            return -1;
        }
        else {
            return hash_entry;
        }
    }*/

    /*kmer_handle_t *init_kmers(char *dataD) {
        kmer_handle_t *handle = malloc(sizeof(kmer_handle_t));

        kmer_memory_image_t *image;

        char file[300];
        strcpy(file,dataD);
        strcat(file,"/function.index");
        handle->function_array = load_functions(file);

        strcpy(file,dataD);
        strcat(file,"/otu.index");
        handle->otu_array      = load_otus(file);

        File fileM = new File(dataD, "kmer.table.mem_map");

            int fd;
            if ((fd = open(fileM, O_RDONLY)) == -1) {
                perror("open");
                exit(1);
            }

            struct stat sbuf;
            if (stat(fileM, &sbuf) == -1) {
                fprintf(stderr, "stat %s failed: %s\n", fileM, strerror(errno));
                exit(1);
            }
            unsigned long long file_size = sbuf.st_size;

            int flags = MAP_SHARED;
            #ifdef MAP_POPULATE
            flags |= MAP_POPULATE;
            #endif

            image = (kmer_memory_image_t *) mmap((caddr_t)0, file_size, PROT_READ, flags, fd, 0);

            if (image == (kmer_memory_image_t *)(-1)) {
                fprintf(stderr, "mmap of kmer_table %s failed: %s\n", fileM, strerror(errno));
                exit(1);
            }

            if (image->version != (long long) VERSION) {
                fprintf(stderr, "Version mismatch for file %s: file has %lld code has %lld\n", 
                        fileM, image->version, (long long) VERSION);
                exit(1);
            }

            if (image->entry_size != (unsigned long long) sizeof(sig_kmer_t)) {
                fprintf(stderr, "Version mismatch for file %s: file has entry size %lld code has %lld\n",
                        fileM, image->entry_size, (unsigned long long) sizeof(sig_kmer_t));
                exit(1);
            }

            size_hash = image->num_sigs;
            handle->num_sigs = size_hash;
            handle->kmer_table = (sig_kmer_t *) (image + 1);

            if (file_size != ((sizeof(sig_kmer_t) * image->num_sigs) + sizeof(kmer_memory_image_t))) {
                fprintf(stderr, "Version mismatch for file %s: file size does not match\n", fileM);
                exit(1);
            }

            fprintf(stderr, "Set size_hash=%lld from file size %lld\n", size_hash, file_size);

        return handle;
    }*/

    /*void advance_past_ambig(unsigned char **p,unsigned char *bound) {

        if (K == 5) {
            while (((*p) < bound) &&
                    ((*(*p) == 20)     || 
                            (*((*p)+1) == 20) || 
                            (*((*p)+2) == 20) || 
                            (*((*p)+3) == 20) || 
                            (*((*p)+4) == 20) )) {
                (*p)++;
            }
        }
        else {   //  ##### ASSUMING K == 8 #### 
            int bad = 1;
            while ((*p < bound) && (bad == 1)) {
                bad = 0;
                if      (*((*p)+7) == 20) {
                    bad = 1;
                    (*p) += 8;
                }
                else if (*((*p)+6) == 20) {
                    bad = 1;
                    (*p) += 7;
                }
                else if (*((*p)+5) == 20) {
                    bad = 1;
                    (*p) += 6;
                }
                else if (*((*p)+4) == 20) {
                    bad = 1;
                    (*p) += 5;
                }
                else if (*((*p)+3) == 20) {
                    bad = 1;
                    (*p) += 4;
                }
                else if (*((*p)+2) == 20) {
                    bad = 1;
                    (*p) += 3;
                }
                else if (*((*p)+1) == 20) {
                    bad = 1;
                    (*p) += 2;
                }
                else if (*((*p)+0) == 20) {
                    bad = 1;
                    (*p) += 1;
                }
            } 
        }
    }

    void display_hits(FILE *fh) {
        fprintf(fh, "hits: ");
        int i;
        for (i=0; (i < num_hits); i++) {
            fprintf(fh, "%d/%f/%d ", hits[i].from0_in_prot,hits[i].function_wt,hits[i].fI);
        }
        fprintf(fh, "\n");
    }


    void process_set_of_hits(kmer_handle_t *kmersH, FILE *fh) {
        int fI_count = 0;
        float weighted_hits = 0;
        int last_hit=0;
        int i=0;
        while (i < num_hits) {
            if (hits[i].fI == current_fI) {
                last_hit = i;
                fI_count++;
                weighted_hits += hits[i].function_wt;
            }
            i++;
        }  if ((fI_count >= min_hits) && (weighted_hits >= min_weighted_hits)) {
            fprintf(fh, "CALL\t%d\t%d\t%d\t%d\t%s\t%f\n",hits[0].from0_in_prot,
                    hits[last_hit].from0_in_prot+(K-1),
                    fI_count,
                    current_fI,
                    kmersH->function_array[current_fI],
                    weighted_hits);

            if (debug > 1) {
                fprintf(fh, "after-call: ");
                display_hits(fh);
            }
            // once we have decided to call a region, we take the kmers for fI and
            //add them to the counts maintained to assign an OTU to the sequence 
            for (i=0; (i <= last_hit); i++) {
                if (hits[i].fI == current_fI) {
                    int j;
                    for (j=0; (j < num_oI) && (oI_counts[j].oI != hits[i].oI); j++) {}
                    if (j == num_oI) {
                        if (num_oI == OI_BUFSZ) {
                            j--;   // we overwrite the last entry
                        }
                        else
                            num_oI++;
                        oI_counts[j].oI    = hits[i].oI;
                        oI_counts[j].count = 1;
                    }
                    else {
                        oI_counts[j].count++;
                    }
                    // now we bubble the count back, allowing it to establish
                    while ((j > 0) && (oI_counts[j-1].count <= oI_counts[j].count)) {
                        int oI_tmp    = oI_counts[j-1].oI;
                        int count_tmp = oI_counts[j-1].count;
                        oI_counts[j-1].oI    = oI_counts[j].oI;
                        oI_counts[j-1].count = oI_counts[j].count;
                        oI_counts[j].oI    = oI_tmp;
                        oI_counts[j].count = count_tmp;
                        j--;
                    }
                }
            }
        }

        if ((hits[num_hits-2].fI != current_fI) && (hits[num_hits-2].fI == hits[num_hits-1].fI)) {
            current_fI = hits[num_hits-1].fI;
            // now copy the last two entries to the start of the hits array.  Sorry this is so clumsy
            hits[0].oI               = hits[num_hits-2].oI;
            hits[0].from0_in_prot    = hits[num_hits-2].from0_in_prot;
            hits[0].avg_off_from_end = hits[num_hits-2].avg_off_from_end;
            hits[0].fI               = hits[num_hits-2].fI;
            hits[0].function_wt      = hits[num_hits-2].function_wt;

            hits[1].oI               = hits[num_hits-1].oI;
            hits[1].from0_in_prot    = hits[num_hits-1].from0_in_prot;
            hits[1].avg_off_from_end = hits[num_hits-1].avg_off_from_end;
            hits[1].fI               = hits[num_hits-1].fI;
            hits[1].function_wt      = hits[num_hits-1].function_wt;
            num_hits                 = 2;
        }
        else {
            num_hits = 0;
        }
    }

    void gather_hits(int ln_DNA, char strand,int prot_off,char *pseq,
            unsigned char *pIseq, kmer_handle_t *kmersH, FILE *fh) {

        if (debug >= 3) {
            fprintf(fh, "translated: %c\t%d\t%s\n",strand,prot_off,pseq);
        }

        unsigned char *p = pIseq;
        // pseq and pIseq are the same length 

        unsigned char *bound = pIseq + strlen(pseq) - K;
        advance_past_ambig(&p,bound);
        unsigned long long encodedK=0;
        if (p < bound) {
            encodedK = encoded_kmer(p);
        }
        while (p < bound) {
            long long  where = lookup_hash_entry(kmersH->kmer_table,encodedK);
            if (where >= 0) {
                sig_kmer_t *kmers_hash_entry = &(kmersH->kmer_table[where]);
                int avg_off_end = kmers_hash_entry->avg_from_end;
                int fI        = kmers_hash_entry->function_index;
                int oI          = kmers_hash_entry->otu_index;
                float f_wt      = kmers_hash_entry->function_wt;
                if (debug >= 1) {
                    fprintf(fh, "HIT\t%ld\t%lld\t%d\t%d\t%0.3f\t%d\n",p-pIseq,encodedK,avg_off_end,fI,f_wt,oI);
                }

                if ((num_hits > 0) && (hits[num_hits-1].from0_in_prot + max_gap) < (p-pIseq)) {
                    if (num_hits >= min_hits) {
                        process_set_of_hits(kmersH, fh);
                    }
                    else {
                        num_hits = 0;
                    }
                }

                if (num_hits == 0) {
                    current_fI = fI;   // if this is the first, set the current_fI
                }

                if ((! order_constraint) || (num_hits == 0) ||
                        ((fI == hits[num_hits-1].fI) &&
                                (abs(((p-pIseq) - hits[num_hits-1].from0_in_prot) - 
                                        (hits[num_hits-1].avg_off_from_end - avg_off_end)
                                        ) <= 20))) {
                    // we have a new hit, so we add it to the global set of hits
                    hits[num_hits].oI = oI;
                    hits[num_hits].fI = fI;
                    hits[num_hits].from0_in_prot = p-pIseq;
                    hits[num_hits].avg_off_from_end = avg_off_end;
                    hits[num_hits].function_wt = f_wt;
                    if (num_hits < MAX_HITS_PER_SEQ - 2) 
                        num_hits++;
                    if (debug > 1) {
                        fprintf(fh, "after-hit: ");
                        display_hits(fh);
                    }
                    if ((num_hits > 1) && (current_fI != fI) &&           // if we have a pair of new fIs, it is time to
                            (hits[num_hits-2].fI == hits[num_hits-1].fI)) {   // process one set and initialize the next
                        process_set_of_hits(kmersH, fh);
                    }
                }
            }
            p++;
            if (p < bound) {
                if (*(p+K-1) < 20) {
                    encodedK = ((encodedK % CORE) * 20L) + *(p+K-1);
                }
                else {
                    p += K;
                    advance_past_ambig(&p,bound);
                    if (p < bound) {
                        encodedK = encoded_kmer(p);
                    }
                }
            }    
        }
        if (num_hits >= min_hits) {
            process_set_of_hits(kmersH, fh);
        }
        num_hits = 0;
    }*/

    /*void tabulate_otu_data_for_contig(FILE *fh) {
        int i;
        fprintf(fh, "OTU-COUNTS\t%s[%d]",current_id,current_length_contig);
        for (i=0; (i < num_oI); i++) {
            fprintf(fh, "\t%d-%d",oI_counts[i].count,oI_counts[i].oI);
        }
        fprintf(fh, "\n");
        num_oI = 0;
    }

    void process_aa_seq(char *id,char *pseq,kmer_handle_t *kmersH, FILE *fh) {
        static unsigned char *pIseq = 0;
        if (pIseq == 0)
        {
            pIseq = malloc(MAX_SEQ_LEN / 3);
        }
        //static unsigned char pIseq[MAX_SEQ_LEN / 3];

        strcpy(current_id,id);
        fprintf(fh, "PROTEIN-ID\t%s\n",id);
        int ln = strlen(pseq);
        current_length_contig = ln;
        current_strand        = '+';
        current_prot_off      = 0;
        int i;
        for (i=0; (i < ln); i++)
            pIseq[i] = to_amino_acid_off(*(pseq+i));
        gather_hits(ln,'+',0,pseq,pIseq,kmersH,fh);  
        tabulate_otu_data_for_contig(fh);
    }

    void process_seq(String id, String data, kmer_handle[] kmersH, PrintWriter pw) {

        int codonCount = data.length() / 3 + 1;
        char[] pseq = new char[codonCount];
        byte[] pIseq = new byte[codonCount];

        String current_id = id;
        int ln = data.length();
        current_length_contig = ln;
        pw.println(String.format("processing %s[%d]",id,ln));
        for (int i = 0; i < 3; i++) {
            translate(data,i,pseq,pIseq);
            current_strand   = '+';
            current_prot_off = i;
            pw.println(String.format("TRANSLATION\t%s\t%d\t%c\t%d",current_id,
                    current_length_contig,
                    current_strand,
                    current_prot_off));
            gather_hits(ln,'+',i,pseq,pIseq,kmersH, pw);
        }
        char[] cdata = rev_comp(data);
        for (int i=0; (i < 3); i++) {
            translate(cdata,i,pseq,pIseq);
            current_strand   = '-';
            current_prot_off = i;
            pw.println(String.format("TRANSLATION\t%s\t%d\t%c\t%d",current_id,
                    current_length_contig,
                    current_strand,
                    current_prot_off));
            gather_hits(ln,'-',i,pseq,pIseq,kmersH, pw);
        }
        tabulate_otu_data_for_contig(pw);
    }

    int main(int argc,char *argv[]) {
        int c;
        char *past;
        char file[300];
        int is_server = 0;
        unsigned short port;

        while ((c = getopt (argc, argv, "ad:s:wD:m:g:OM:l:")) != -1) {
            switch (c) {
            case 'a':
                aa = 1;
                break;
            case 'd':
                debug = strtol(optarg,&past,0);
                break;
            case 'l':
                port = atoi(optarg);
                is_server = 1;
                break;
            case 'm':
                min_hits = strtol(optarg,&past,0);
                break;
            case 'M':
                min_weighted_hits = strtol(optarg,&past,0);
                break;
            case 'O':
                order_constraint = 1;
                break;
            case 'g':
                max_gap = strtol(optarg,&past,0);
                break;
            case 'D':
                strcpy(file,optarg);
                break;
            case 's':
                size_hash = strtol(optarg,&past,0);
                break;
            case 'w':
                write_mem_map = 1;
                break;
            default:
                fprintf(stderr,"arguments: [-a] [-d level] [-s hash-size] [-w] [-m min_hits] -D DataDir \n");
                abort ();
            }
        }

        kmer_handle_t *kmersH = init_kmers(file);

        run_from_filehandle(kmersH, stdin, stdout);
        return 0;
    }

    void run_from_filehandle(kmer_handle_t *kmersH, FILE *fh_in, FILE *fh_out)
    {
        static char *data = 0;
        if (data == 0)
        {
            data = malloc(MAX_SEQ_LEN);
        }

        // static char data[MAX_SEQ_LEN];
        int got_gt = 0;
        int i;
        char *p;
        char id[2000];

        while (((!got_gt && (fscanf(fh_in,">%s",id) == 1)) ||
                (got_gt && (fscanf(fh_in,"%s",id) == 1)))) {
            while (getc(fh_in) != '\n')
                ;
            if ((id[0] == 'F') &&
                    (id[1] == 'L') &&
                    (id[2] == 'U') &&
                    (id[3] == 'S') &&
                    (id[4] == 'H')) {   // ugly compare, since \n is included in id 

                fprintf(fh_out, "//\n");
                got_gt = 0;
            }
            else {

                for (p=data; ((i = getc(fh_in)) != -1) && (i != '>');) {
                    if ((i != ' ') && (i != '\n'))
                        *(p++) = toupper(i);
                }

                if (i == '>')
                    got_gt = 1;
                else
                    got_gt = 0;

                *p=0;
                if ((p-data) > MAX_SEQ_LEN) {
                    fprintf(stderr,"The contig size exceeds %d; bump MAX_SEQ_LEN\n",MAX_SEQ_LEN);
                    exit(1);
                }

                if (! aa)
                    process_seq(id,data,kmersH, fh_out);
                else
                    process_aa_seq(id,data,kmersH, fh_out);
                fflush(fh_out);
            }
        }

        if (debug >= 2)
            fprintf(fh_out, "tot_lookups=%d retry=%d\n",tot_lookups,retry);
    }*/

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("Usage: <program> <kmer-table> <contigs-fasta>");
            System.exit(1);
        }
        FastaReader fr = new FastaReader(new File(args[1]));
        Map<Long, QueryKmer> queryMap = new HashMap<Long, QueryKmer>();
        while (true) {
            String[] entry = fr.read();
            if (entry == null)
                break;
            char[] seq = entry[1].toCharArray();
            int len = seq.length / 3 + 1;
            char[] pseq = new char[len];
            byte[] pIseq = new byte[len];
            translate(seq, 0, pseq, pIseq);
            for (int i = 0; i < len - 10; i++) {
                long value = encoded_kmer(pIseq, i);
                if (value < 0)
                    continue;
                QueryKmer qk = queryMap.get(value);
                if (qk == null) {
                    qk = new QueryKmer();
                    qk.value = value;
                    qk.posList = new ArrayList<QueryPos>(1);
                    queryMap.put(value, qk);
                }
                QueryPos qp = new QueryPos();
                qp.queryId = entry[0];
                qp.pos = i;
                qk.posList.add(qp);
            }
        }
        List<QueryKmer> values = new ArrayList<QueryKmer>(queryMap.values());
        System.out.println("Value count: " + values.size());
        InputStream is = new BufferedInputStream(new FileInputStream(new File(args[0])));
        long t1 = System.currentTimeMillis();
        int kmersFound = 0;
        try {
            long numSigs = readLongLE(is);
            long entrySize = readLongLE(is);
            long version = readLongLE(is);
            System.out.println("numSigs=" + numSigs + ", entrySize=" + entrySize + ", version=" + version);
            // Update hash-codes in queries and sort them by these hash-codes
            for (QueryKmer qk : values) {
                qk.hashCode = qk.value % numSigs;
            }
            Collections.sort(values, new Comparator<QueryKmer>() {
                @Override
                public int compare(QueryKmer o1, QueryKmer o2) {
                    int ret = Long.compare(o1.hashCode, o2.hashCode);
                    if (ret == 0) {
                        ret = Long.compare(o1.value, o2.value);
                    }
                    return ret;
                }
            });
            // Now we go along hash table and along query hash-codes
            long curHashCode = 0;
            int curQueryPos = 0;
            Map<Long, QueryKmer> inProgress = new HashMap<Long, QueryKmer>();
            int fraction = 0;
            while (curQueryPos < values.size() || inProgress.size() > 0) {
                if (inProgress.size() == 0) {
                    QueryKmer qk = values.get(curQueryPos);
                    long neededHashCode = qk.hashCode;
                    // Let's push all queries with the same 
                    inProgress.put(qk.value, qk);
                    curQueryPos++;
                    while (curQueryPos < values.size()) {
                        qk = values.get(curQueryPos);
                        if (qk.hashCode != neededHashCode) {
                            break;
                        }
                        inProgress.put(qk.value, qk);
                        curQueryPos++;
                    }
                    long bytesToSkip = entrySize * (long)(neededHashCode - curHashCode);
                    if (bytesToSkip > 0) {
                        long bytesLeft = bytesToSkip;
                        while (bytesLeft > 0) {
                            long ret = is.skip(bytesLeft);
                            if (ret < 0)
                                break;
                            bytesLeft -= ret;
                        }
                        if (bytesLeft != 0) {
                            throw new IllegalStateException("Error skipping " + bytesToSkip + " bytes");
                        }
                    }
                    curHashCode = neededHashCode;
                }
                long whichKmer = readLongLE(is);
                int otuIndex = readIntLE(is);
                int avgFromEnd = readIntLE(is);
                int functionIndex = readIntLE(is);
                float functionWt = readFloatLE(is);
                long slot = whichKmer % numSigs;
                if (whichKmer > MAX_ENCODED) {
                    inProgress.clear();
                } else {
                    if (inProgress.containsKey(whichKmer)) {
                        QueryKmer qk = inProgress.remove(whichKmer);
                        System.out.println("[" + curHashCode + "] whichKmer=" + whichKmer + " (" + slot + "), otuIndex=" + otuIndex + ", " +
                                "avgFromEnd=" + avgFromEnd + ", functionIndex=" + functionIndex + ", " +
                                "functionWt=" + functionWt);
                        sig_kmer hit = new sig_kmer();
                        hit.which_kmer = whichKmer;
                        hit.otu_index = otuIndex;
                        hit.avg_from_end = avgFromEnd;
                        hit.function_index = functionIndex;
                        hit.function_wt = functionWt;
                        qk.hit = hit;
                        kmersFound++;
                    }
                }
                curHashCode++;
                int newFraction = (int)(1000.0 * ((double)curHashCode / (double)numSigs));
                if (newFraction != fraction) {
                    fraction = newFraction;
                    System.out.println("Processed: " + (fraction / 10.0) + "%, time=" +
                            (System.currentTimeMillis() - t1) + " ms.");
                }
            }
            System.out.println("curQueryPos=" + curQueryPos + " (size=" + values.size() + "), n=" + curHashCode);
        } finally {
            is.close();
            System.out.println("Kmers found: " + kmersFound);
            System.out.println("Time: " + (System.currentTimeMillis() - t1) + " ms.");
        }
    }

    public static int readIntLE(InputStream is) throws IOException {
        int ch1 = is.read();
        int ch2 = is.read();
        int ch3 = is.read();
        int ch4 = is.read();
        if ((ch1 | ch2 | ch3 | ch4) < 0)
            throw new EOFException();
        return ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
    }

    public static long readLongLE(InputStream is) throws IOException {
        return (((long)(is.read() & 255) << 0) +
                ((long)(is.read() & 255) << 8) +
                ((long)(is.read() & 255) << 16) +
                ((long)(is.read() & 255) << 24) +
                ((long)(is.read() & 255) << 32) +
                ((long)(is.read() & 255) << 40) +
                ((long)(is.read() & 255) << 48) +
                ((long)is.read() << 56));
    }

    public static float readFloatLE(InputStream is) throws IOException {
        return Float.intBitsToFloat(readIntLE(is));
    }

    public static class FastaReader {
        BufferedReader br;
        String str1 = null;

        public FastaReader(File f) {
            try {
                br = new BufferedReader(new FileReader(f));
            } catch(Exception ex) {
                throw new RuntimeException("Wrong file name: " + f, ex);
            }
        }

        public FastaReader(Reader r) {
            br = new BufferedReader(r);
        }

        public static Map<String, String> readFromFile(File f) {
            FastaReader fr = new FastaReader(f);
            Map<String, String> ret = fr.readAll();
            fr.close();
            return ret;
        }

        public Map<String, String> readAll() {
            Map<String, String> ret = new LinkedHashMap<String, String>();
            while (true) {
                String[] entry = read();
                if (entry == null)
                    break;
                ret.put(entry[0], entry[1]);
            }
            return ret;
        }

        public String[] read() {
            if(br == null)
                return null;
            String protName = null;
            String protSeq = null;
            String protDescr = null;
            try {
                if(str1 == null)
                    str1 = br.readLine();
                for(;;) {
                    if(str1 == null)
                        return null;
                    String str2 = str1.trim();
                    if(str2.length() > 1) {
                        if (str2.charAt(0) == '>' && str2.substring(1).trim().length() > 0) {
                            StringTokenizer st = new StringTokenizer(str2.substring(1), " \t");
                            protName = st.nextToken();
                            StringBuilder descr = new StringBuilder();
                            while (st.hasMoreTokens()) {
                                if (descr.length() > 0)
                                    descr.append(" ");
                                descr.append(st.nextToken());
                            }
                            protDescr = descr.toString();
                            break;
                        } else {
                            throw new IllegalStateException("Wrong caption line: " + str2);
                        }
                    }
                    str1 = br.readLine();
                }
                StringBuilder sb = new StringBuilder();
                for(;;) {
                    str1 = br.readLine();
                    if(str1 == null || str1.trim().startsWith(">")) {
                        throw new IllegalStateException("No sequence for caption: " + protName);
                    }
                    if(str1.trim().length() > 0)
                        break;
                }
                for(;;) {
                    sb.append(str1);
                    str1 = br.readLine();
                    if(str1 == null || str1.trim().startsWith(">"))
                        break;
                }
                protSeq = sb.toString();
                if (protSeq.length()==0) {
                    throw new IllegalStateException("No sequence for caption: " + protName);
                }
                return new String[] {protName, protSeq, protDescr};
            } catch (IllegalStateException ex) {
                throw ex;
            } catch (Exception ex) {
                throw new IllegalStateException(ex);
            }
        }
        public void close() {
            if(br == null)
                return;
            try {
                br.close();
            }catch(Exception ex) {
                System.err.println("WARNING: couldn't close fasta reader, ignored");
            }
            br = null;
        }
    }
    
    public static class QueryKmer {
        long value;
        long hashCode;
        List<QueryPos> posList;
        sig_kmer hit;
    }
    
    public static class QueryPos {
        String queryId;
        int pos;
    }
}