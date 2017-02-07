package kmergutsjava;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

import us.kbase.common.service.UObject;

import com.fasterxml.jackson.core.type.TypeReference;

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

    public static final int K = 8;
    public static final int MAX_SEQ_LEN = 500000000;

    public static final long CORE = 20L*20L*20L*20L*20L*20L*20L;

    public static final long MAX_ENCODED = CORE*20L; 

    public static final char[] GENETIC_CODE = {
            'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
            'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
            'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
            '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
    };

    public static final char[] PROT_ALPHA = {
            'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
    };

    public static final int VERSION = 1;

    public static final int MAX_HITS_PER_SEQ = 40000;

    public static final int OI_BUFSZ = 5;
    
    boolean aa = false;
    boolean orderConstraint = false;
    int minHits = 5;
    int minWeightedHits = 0;
    int maxGap  = 200;

    public static final int MAX_FUNC_OI_INDEX = 1000000;
    public static final int MAX_FUNC_OI_VALS  = 100000000;

    /* =========================== end of reduction global variables ================= */

    public static byte toAminoAcidOff(char c) {
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


    public static char[] revComp(char[] data) {
        int n = data.length;
        char[] cdata = new char[n];
        int p  = n - 1;
        int pc = 0;
        while (n-- > 0) {
            cdata[pc++] = compl(data[p--]);
        }
        return cdata;
    }

    public static long encodedKmer(byte[] data, int pos) {
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

    public static int dnaChar(char c)
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
            int c1 = dnaChar(seq[i++]);
            int c2 = dnaChar(seq[i++]);
            int c3 = dnaChar(seq[i++]);
            if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
                int I = (c1 * 16) + (c2 * 4) + c3;
                char protC = GENETIC_CODE[I];
                pseq[p] = protC;
                pIseq[p] = toAminoAcidOff(protC);
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

    private static List<String> loadIndexedArray(File file) throws Exception {
        BufferedReader br;
        if (file.getName().endsWith(".gz")) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                    new BufferedInputStream(new FileInputStream(file)))));
        } else {
            br = new BufferedReader(new FileReader(file));
        }

        List<String> ret = new ArrayList<String>();
        for (int linePos = 0; ; linePos++) {
            String line = br.readLine();
            if (line == null)
                break;
            int tabPos = line.indexOf('\t');
            int index = Integer.parseInt(line.substring(0, tabPos));
            if (linePos != index) {
                throw new IllegalStateException("Your index must be dense and in order (see line " + linePos + ")");
            }
            ret.add(line.substring(tabPos + 1));
        }
        br.close();
        return ret;
    }

    private List<String> loadFunctions(File file) throws Exception {
        return loadIndexedArray(file);
    }

    /*List<String> load_otus(File file) {
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
    }*/


    public int processSetOfHits(List<Hit> hits, List<String> functionArray, 
            int currentFI, List<OtuCount> oICounts, PrintWriter out) {
        int fICount = 0;
        float weightedHits = 0;
        int lastHit = 0;
        for (int i = 0; i < hits.size(); i++) {
            if (hits.get(i).fI == currentFI) {
                lastHit = i;
                fICount++;
                weightedHits += hits.get(i).functionWt;
            }
            i++;
        }
        if ((fICount >= minHits) && (weightedHits >= minWeightedHits)) {
            out.println(String.format("CALL\t%d\t%d\t%d\t%d\t%s\t%f",
                    hits.get(0).from0InProt,
                    hits.get(lastHit).from0InProt + (K-1),
                    fICount,
                    currentFI,
                    functionArray.get(currentFI),
                    weightedHits));

            // once we have decided to call a region, we take the kmers for fI and
            //add them to the counts maintained to assign an OTU to the sequence
            for (int i=0; i <= lastHit; i++) {
                if (hits.get(i).fI == currentFI) {
                    int j;
                    for (j=0; (j < oICounts.size()) && (oICounts.get(j).oI != hits.get(i).oI); j++) {}
                    if (j == oICounts.size()) {
                        if (oICounts.size() == OI_BUFSZ) {
                            j--;   // we overwrite the last entry
                        }
                        else {
                            oICounts.add(new OtuCount());
                        }
                        oICounts.get(j).oI = hits.get(i).oI;
                        oICounts.get(j).count = 1;
                    }
                    else {
                        oICounts.get(j).count++;
                    }
                    // now we bubble the count back, allowing it to establish
                    while ((j > 0) && (oICounts.get(j-1).count <= oICounts.get(j).count)) {
                        OtuCount tmp = oICounts.get(j-1);
                        oICounts.set(j-1, oICounts.get(j));
                        oICounts.set(j, tmp);
                        j--;
                    }
                }
            }
        }
        int numHits = hits.size();
        if ((hits.get(numHits-2).fI != currentFI) && (hits.get(numHits - 2).fI == hits.get(numHits - 1).fI)) {
            currentFI = hits.get(numHits - 1).fI;
            // now copy the last two entries to the start of the hits array.  Sorry this is so clumsy
            hits.set(0, hits.get(numHits - 2));
            hits.set(1, hits.get(numHits - 1));
            for (int i = numHits - 1; i >= 2; i--) {
                hits.remove(i);
            }
        }
        else {
            hits.clear();
        }
        return currentFI;
    }

    public void gatherHits(int ln_DNA, char strand, int frame, 
            List<Hit> allHits, List<String> functionArray, 
            List<OtuCount> oICounts, PrintWriter pw) {
        Collections.sort(allHits, new Comparator<Hit>() {
            @Override
            public int compare(Hit o1, Hit o2) {
                return Integer.compare(o1.from0InProt, o2.from0InProt);
            }
        });
        List<Hit> hits = new ArrayList<Hit>();
        int currentFI = 0;
        for (Hit ph : allHits) {
                int avgOffEnd = ph.avgOffFromEnd;
                int fI = ph.fI;

                if ((hits.size() > 0) && (hits.get(hits.size()-1).from0InProt + maxGap) < ph.from0InProt) {
                    if (hits.size() >= minHits) {
                        currentFI = processSetOfHits(hits, functionArray, currentFI, oICounts, pw);
                    }
                    else {
                        hits.clear();
                    }
                }

                if (hits.isEmpty()) {
                    currentFI = fI;   // if this is the first, set the current_fI
                }

                if ((!orderConstraint) || (hits.size() == 0) ||
                        ((fI == hits.get(hits.size()-1).fI) &&
                                (Math.abs(((ph.from0InProt) - hits.get(hits.size()-1).from0InProt) - 
                                        (hits.get(hits.size()-1).avgOffFromEnd - avgOffEnd)
                                        ) <= 20))) {
                    // we have a new hit, so we add it to the global set of hits
                    if (hits.size() < MAX_HITS_PER_SEQ - 2) {
                        hits.add(ph);
                    }
                    if ((hits.size() > 1) && (currentFI != fI) &&           // if we have a pair of new fIs, it is time to
                            (hits.get(hits.size()-2).fI == hits.get(hits.size()-1).fI)) {   // process one set and initialize the next
                        currentFI = processSetOfHits(hits, functionArray, currentFI, oICounts, pw);
                    }
                }
        }
        if (hits.size() >= minHits) {
            processSetOfHits(hits, functionArray, currentFI, oICounts, pw);
        }
    }

    void tabulateOtuDataForContig(String currentId, int currentLengthContig, 
            List<OtuCount> oICounts, PrintWriter pw) {
        pw.print(String.format("OTU-COUNTS\t%s[%d]", currentId, currentLengthContig));
        for (OtuCount oICount : oICounts) {
            pw.print(String.format("\t%d-%d", oICount.count, oICount.oI));
        }
        pw.println();
        oICounts.clear();
    }

    public void processAASeq(String id, int proteinLen, Map<HitContainerKey, HitContainer> hitCnts,
            List<String> functionArray, PrintWriter pw) {
        List<OtuCount> oICounts = new ArrayList<OtuCount>();
        pw.println(String.format("PROTEIN-ID\t%s", id));
        HitContainerKey key = new HitContainerKey();
        key.queryId = id;
        key.strand = '+';
        key.frame = 0;
        gatherHits(proteinLen, '+', 0, hitCnts.get(key).hits, functionArray, oICounts, pw);  
        tabulateOtuDataForContig(id, proteinLen, oICounts, pw);
    }

    void processSeq(String id, int contigLen, Map<HitContainerKey, HitContainer> hitCnts, 
            List<String> functionArray, PrintWriter pw) {
        List<OtuCount> oICounts = new ArrayList<OtuCount>();
        pw.println(String.format("processing %s[%d]", id, contigLen));
        char[] strands = {'+', '-'};
        for (char strand : strands) {
            for (int frame = 0; frame < 3; frame++) {
                pw.println(String.format("TRANSLATION\t%s\t%d\t%c\t%d", id,
                        contigLen,
                        strand,
                        frame));
                HitContainerKey key = new HitContainerKey();
                key.queryId = id;
                key.strand = strand;
                key.frame = frame;
                gatherHits(contigLen, strand, frame, hitCnts.get(key).hits, functionArray, 
                        oICounts, pw);
            }
        }
        tabulateOtuDataForContig(id, contigLen, oICounts, pw);
    }

    /*int main(int argc,char *argv[]) {
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
                order_constraint = true;
                break;
            case 'g':
                max_gap = strtol(optarg,&past,0);
                break;
            case 'D':
                strcpy(file,optarg);
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
        try {
            String kmerTableDir = args[0];
            String contigsPath = args[1];
            KmerGutsJava instance = new KmerGutsJava();
            PrintWriter pw = new PrintWriter(System.out);
            instance.run(new File(kmerTableDir), new File(contigsPath), pw);
            pw.flush();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    public void run(File kmerTableDir, File contigsFile, PrintWriter pw) throws Exception {
        File kmerTableFile = new File(kmerTableDir, "kmer.table.mem_map");
        File kmerTableGzFile = new File(kmerTableDir, kmerTableFile.getName() + ".gz");
        if (kmerTableGzFile.exists()) {
            kmerTableFile = kmerTableGzFile;
        }
        File hitsFile = new File("/kb/module/work/hits.json");
        //File hitsFile = new File("test_local/workdir/hits.json");
        {
            List<HitContainer> hits = lookup(kmerTableFile, contigsFile);
            System.out.println("Writing hits to JSON file");
            UObject.getMapper().writeValue(hitsFile, hits);
            System.out.println("Done");
        }
        List<HitContainer> hits = UObject.getMapper().readValue(hitsFile, new TypeReference<List<HitContainer>>() {});
        Map<HitContainerKey, HitContainer> hitCnts = 
                new LinkedHashMap<HitContainerKey, HitContainer>();
        for (HitContainer cnt : hits) {
            hitCnts.put(cnt.key, cnt);
        }
        File functionIndexFile = new File(kmerTableDir, "function.index");
        File functionIndexGzFile = new File(kmerTableDir, functionIndexFile.getName() + ".gz");
        if (functionIndexGzFile.exists()) {
            functionIndexFile = functionIndexGzFile;
        }
        List<String> functionArray = loadFunctions(functionIndexFile);
        FastaReader fr = new FastaReader(contigsFile);
        while (true) {
            String[] entry = fr.read();
            if (entry == null)
                break;
            String id = entry[0];
            int seqLen = entry[1].length();
            if (aa) {
                processAASeq(id, seqLen, hitCnts, functionArray, pw);
            } else {
                processSeq(id, seqLen, hitCnts, functionArray, pw);
            }
            pw.flush();
        }
    }
    
    private static void addKmers(String id, char strand, int frame, char[] pseq,
            byte[] pIseq, List<QueryKmer> queryList, List<HitContainer> hitCnts) {
        HitContainerKey hcKey = new HitContainerKey();
        hcKey.queryId = id;
        hcKey.strand = strand;
        hcKey.frame = frame;
        HitContainer hitCnt = new HitContainer();
        hitCnt.key = hcKey;
        hitCnt.hits = new ArrayList<Hit>();
        hitCnts.add(hitCnt);
        for (int i = 0; i < pIseq.length - K; i++) {
            long value = encodedKmer(pIseq, i);
            if (value < 0)
                continue;
            QueryKmer qk = new QueryKmer();
            qk.value = value;
            qk.protPos = i;
            qk.hitCnt = hitCnt;
            queryList.add(qk);
        }
    }
    
    private List<HitContainer> lookup(File kmerTableFile, File contigsPath) throws Exception {
        List<QueryKmer> values = new ArrayList<QueryKmer>();
        List<HitContainer> hitCnts = new ArrayList<HitContainer>();
        FastaReader fr = new FastaReader(contigsPath);
        while (true) {
            String[] entry = fr.read();
            if (entry == null)
                break;
            String id = entry[0];
            char[] seq = entry[1].toCharArray();
            if (aa) {
                byte[] pIseq = new byte[seq.length];
                for (int i = 0; i < seq.length; i++) {
                    pIseq[i] = toAminoAcidOff(seq[i]);
                }
            } else {
                int len = seq.length / 3 + 1;
                char[] pseq = new char[len];
                byte[] pIseq = new byte[len];
                for (int frame = 0; frame < 3; frame++) {
                    translate(seq, frame, pseq, pIseq);
                    addKmers(id, '+', frame, pseq, pIseq, values, hitCnts);
                }
                char[] complSeq = revComp(seq);
                for (int frame = 0; frame < 3; frame++) {
                    translate(complSeq, frame, pseq, pIseq);
                    addKmers(id, '-', frame, pseq, pIseq, values, hitCnts);
                }
            }
        }
        fr.close();
        System.out.println("Value count: " + values.size());
        InputStream is;
        if (kmerTableFile.getName().endsWith(".gz")) {
            is = new BufferedInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(kmerTableFile))));
        } else {
            is = new BufferedInputStream(new FileInputStream(kmerTableFile));
        }
        long numSigs = readLongLE(is);
        long entrySize = readLongLE(is);
        long version = readLongLE(is);
        System.out.println("numSigs=" + numSigs + ", entrySize=" + entrySize + ", version=" + version);
        // Update hash-codes in queries and sort them by these hash-codes
        for (QueryKmer qk : values) {
            qk.hashCode = qk.value % numSigs;
        }
        // Sort query kmers by hash-codes (to be able to process them in parrallel with
        // processing kmer-table file (sorted joining).
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
            
        long t1 = System.currentTimeMillis();
        int kmersFound = 0;
        int posCount = 0;
        try {
            // Now we go along hash table and along query hash-codes
            long curHashCode = 0;
            int curQueryPos = 0;
            Map<Long, List<QueryKmer>> inProgress = new HashMap<Long, List<QueryKmer>>();
            int fraction = 0;
            while (curQueryPos < values.size() || inProgress.size() > 0) {
                long neededHashCode = curHashCode;
                if (inProgress.size() == 0) {
                    // Update next hash-code
                    QueryKmer qk = values.get(curQueryPos);
                    neededHashCode = qk.hashCode;
                    List<QueryKmer> list = new ArrayList<QueryKmer>(5);
                    list.add(qk);
                    inProgress.put(qk.value, list);
                    curQueryPos++;
                }
                // Let's push all queries with necessary hash-code into progress state
                while (curQueryPos < values.size()) {
                    QueryKmer qk = values.get(curQueryPos);
                    if (qk.hashCode != neededHashCode) {
                        break;
                    }
                    if (inProgress.containsKey(qk.value)) {
                        inProgress.get(qk.value).add(qk);
                    } else {
                        List<QueryKmer> list = new ArrayList<QueryKmer>(3);
                        list.add(qk);
                        inProgress.put(qk.value, list);
                    }
                    curQueryPos++;
                }
                // Let's position kmer-table stream to necessary hash-code
                if (neededHashCode > curHashCode) {
                    long bytesToSkip = entrySize * (long)(neededHashCode - curHashCode);
                    if (bytesToSkip > 0) {
                        long bytesLeft = bytesToSkip;
                        for (int iter = 0; bytesLeft > 0 && iter < 100; iter++) {
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
                if (whichKmer > MAX_ENCODED) {
                    inProgress.clear();
                } else {
                    // Matching kmers from query and from kmer-table
                    if (inProgress.containsKey(whichKmer)) {
                        kmersFound++;
                        //long slot = whichKmer % numSigs;
                        //System.out.println("[" + curHashCode + "] whichKmer=" + whichKmer + " (" + slot + "), otuIndex=" + otuIndex + ", " +
                        //        "avgFromEnd=" + avgFromEnd + ", functionIndex=" + functionIndex + ", " +
                        //        "functionWt=" + functionWt);
                        for (QueryKmer qk : inProgress.remove(whichKmer)) {
                            Hit hit = new Hit();
                            hit.from0InProt = qk.protPos;
                            hit.oI = otuIndex;
                            hit.avgOffFromEnd = avgFromEnd;
                            hit.fI = functionIndex;
                            hit.functionWt = functionWt;
                            qk.hitCnt.hits.add(hit);
                            posCount++;
                        }
                    }
                }
                curHashCode++;
                int newFraction = (int)(100.0 * ((double)curHashCode / (double)numSigs));
                if (newFraction != fraction) {
                    fraction = newFraction;
                    System.out.println("Processed: " + (fraction) + "%, time=" +
                            (System.currentTimeMillis() - t1) + " ms., found-so-far=" + 
                            kmersFound);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            is.close();
        }
        System.out.println("Kmers found: " + kmersFound + " (pos-count=" + posCount + ")");
        System.out.println("Time: " + (System.currentTimeMillis() - t1) + " ms.");
        return hitCnts;
    }

    /*public static void main(String[] args) throws Exception {
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
        RandomAccessFile is = new RandomAccessFile(new File(args[0]), "r");
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
            int fraction = 0;
            for (QueryKmer qk : values) {
                int newFraction = (int)(10000.0 * ((double)qk.hashCode / (double)numSigs));
                if (newFraction != fraction) {
                    fraction = newFraction;
                    System.out.println("Processed: " + (fraction / 100.0) + "%, time=" +
                            (System.currentTimeMillis() - t1) + " ms.");
                }
                long pos = 8 * 3 + qk.hashCode * entrySize;
                is.seek(pos);
                for (long curHashCode = qk.hashCode; ; curHashCode++) {
                    long whichKmer = readLongLE(is);
                    int otuIndex = readIntLE(is);
                    int avgFromEnd = readIntLE(is);
                    int functionIndex = readIntLE(is);
                    float functionWt = readFloatLE(is);
                    long slot = whichKmer % numSigs;
                    if (whichKmer > MAX_ENCODED) {
                        break;
                    } else {
                        if (whichKmer == qk.value) {
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
                            break;
                        }
                    }
                }
            }
        } finally {
            is.close();
            System.out.println("Kmers found: " + kmersFound);
            System.out.println("Time: " + (System.currentTimeMillis() - t1) + " ms.");
        }
    }*/

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
        int ch1 = is.read();
        int ch2 = is.read();
        int ch3 = is.read();
        int ch4 = is.read();
        int ch5 = is.read();
        int ch6 = is.read();
        int ch7 = is.read();
        int ch8 = is.read();
        if ((ch1 | ch2 | ch3 | ch4 | ch5 | ch6 | ch7 | ch8) < 0)
            throw new EOFException();
        return (((long)(ch1 & 255) << 0) +
                ((long)(ch2 & 255) << 8) +
                ((long)(ch3 & 255) << 16) +
                ((long)(ch4 & 255) << 24) +
                ((long)(ch5 & 255) << 32) +
                ((long)(ch6 & 255) << 40) +
                ((long)(ch7 & 255) << 48) +
                ((long)ch8 << 56));
    }

    public static float readFloatLE(InputStream is) throws IOException {
        return Float.intBitsToFloat(readIntLE(is));
    }

    public static int readIntLE(RandomAccessFile is) throws IOException {
        int ch1 = is.read();
        int ch2 = is.read();
        int ch3 = is.read();
        int ch4 = is.read();
        if ((ch1 | ch2 | ch3 | ch4) < 0)
            throw new EOFException();
        return ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
    }

    public static long readLongLE(RandomAccessFile is) throws IOException {
        return (((long)(is.read() & 255) << 0) +
                ((long)(is.read() & 255) << 8) +
                ((long)(is.read() & 255) << 16) +
                ((long)(is.read() & 255) << 24) +
                ((long)(is.read() & 255) << 32) +
                ((long)(is.read() & 255) << 40) +
                ((long)(is.read() & 255) << 48) +
                ((long)is.read() << 56));
    }

    public static float readFloatLE(RandomAccessFile is) throws IOException {
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
        public long value;
        public long hashCode;
        public int protPos;
        public HitContainer hitCnt;
    }
    
    public static class QueryPos {
        public String queryId;
        public char strand;
        public byte frame;
        public int pos;
    }
    
    public static class Hit {
        public int oI;
        public int from0InProt;      /* offset from start of protein sequence */
        public int avgOffFromEnd;  /* average offset from the end */
        public int fI;
        public float functionWt;
    }

    public static class OtuCount {
        int oI;
        int count;
    }
    
    public static class HitContainerKey {
        public String queryId;
        public char strand;
        public int frame;
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + frame;
            result = prime * result
                    + ((queryId == null) ? 0 : queryId.hashCode());
            result = prime * result + strand;
            return result;
        }
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            HitContainerKey other = (HitContainerKey) obj;
            if (frame != other.frame)
                return false;
            if (queryId == null) {
                if (other.queryId != null)
                    return false;
            } else if (!queryId.equals(other.queryId))
                return false;
            if (strand != other.strand)
                return false;
            return true;
        }
    }
    
    public static class HitContainer {
        public HitContainerKey key;
        public List<Hit> hits;
    }
}