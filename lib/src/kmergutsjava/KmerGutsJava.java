package kmergutsjava;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

/*
This is Java version of kmer_guts.c program (http://biocvs.mcs.anl.gov/viewcvs.cgi/Kmers2/kmer_guts.c?revision=1.34&content-type=text%2Fplain)
Usage example:

   kmer_guts -D KmerData < input.contigs

It can also be used to call translated PEGs:

   kmer_guts -D KmerData -a  < amino_acid_sequences.fasta > hits

where KmersData is a directory that must contain

       function.index     [a file of [index,function] pairs]
       kmer.table.mem_map [a binary file with KMER hash table]
----------------

Conceptually, the data associated with each signature Kmer is

        1. the protein kmer
        2. the average offset from the end of the protein
        3. a set of numeric values that include

                function index
                OTU index

The program takes as input a fasta file.
In effect, the program processes a sequence of requests.  The output for each request
is a piece of stdout terminated by a line

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

This code uses a table indicating which K-mers are signatures.
 */
public class KmerGutsJava {
    // Constants
    public static final int K = 8;
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

    /* parameters to main -- accessed globally */
    private boolean aa = false;
    private boolean orderConstraint = false;
    private int minHits = 5;
    private int minWeightedHits = 0;
    private int maxGap  = 200;
    private boolean debug = false;
    private long inputSizeLimit = 10 * 1024 * 1024; // 10 megabases

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

    private static void translate(char[] seq, int off, char[] pseq, byte[] pIseq) {
        int max = seq.length - 3;
        int p = 0;
        for (int i = off; i <= max; ) {
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

    void displayHits(List<Hit> hits, PrintWriter pw) {
        pw.print("hits: ");
        int i;
        for (i=0; i < hits.size(); i++) {
            Hit h = hits.get(i);
            pw.print(String.format("%d/%f/%d ", h.from0InProt, h.functionWt, h.fI));
        }
        pw.println();
    }

    public int processSetOfHits(List<Hit> hits, List<String> functionArray, 
            int currentFI, List<OtuCount> oICounts, PrintWriter pw) {
        int fICount = 0;
        float weightedHits = 0;
        int lastHit = 0;
        for (int i = 0; i < hits.size(); i++) {
            if (hits.get(i).fI == currentFI) {
                lastHit = i;
                fICount++;
                weightedHits += hits.get(i).functionWt;
            }
        }
        if ((fICount >= minHits) && (weightedHits >= minWeightedHits)) {
            pw.println(String.format("CALL\t%d\t%d\t%d\t%d\t%s\t%f",
                    hits.get(0).from0InProt,
                    hits.get(lastHit).from0InProt + (K-1),
                    fICount,
                    currentFI,
                    functionArray.get(currentFI),
                    weightedHits));

            if (debug) {
                pw.print("after-call: ");
                displayHits(hits, pw);
            }

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

                if (debug) {
                    pw.println(String.format("HIT\t%d\t%d\t%d\t%d\t%1.3f\t%d", 
                            ph.from0InProt, 0, avgOffEnd, fI, ph.functionWt, ph.oI));
                }

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
                        if (debug) {
                            pw.print("after-hit: ");
                            displayHits(hits, pw);
                        }
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

    public static void main(String[] args) throws Exception {
        String kmerTableDir = null;
        String queryFastaFile = null;
        String outputFile = null;
        KmerGutsJava instance = new KmerGutsJava();
        try {
            Queue<String> params = new LinkedList<String>(Arrays.asList(args));
            while (!params.isEmpty()) {
                String param = params.poll();
                if(!param.startsWith("-")) {
                    throw new IllegalStateException("Parameter name should start from '-': " + param);
                }
                param = param.substring(1);
                if (param.length() != 1) {
                    throw new IllegalStateException("Unknown parameter: -" + param);
                }
                switch (param.charAt(0)) {
                case 'a':
                    instance.aa = true;
                    break;
                case 'd':
                    instance.debug = true;
                    break;
                case 'm':
                    instance.minHits = Integer.parseInt(params.poll());
                    break;
                case 'M':
                    instance.minWeightedHits = Integer.parseInt(params.poll());
                    break;
                case 'O':
                    instance.orderConstraint = true;
                    break;
                case 'g':
                    instance.maxGap = Integer.parseInt(params.poll());
                    break;
                case 'D':
                    kmerTableDir = params.poll();
                    break;
                case 'q':
                    queryFastaFile = params.poll();
                    break;
                case 'o':
                    outputFile = params.poll();
                    break;
                default:
                    throw new IllegalStateException("Unknown parameter: -" + param);
                }
            }
            if (kmerTableDir == null) {
                throw new IllegalStateException("-D parameter is required");
            }
        } catch (Exception ex) {
            System.out.println("Error: " + ex.getMessage());
            System.out.println("Usage: kmer_guts [options] -D DataDir");
            System.out.println("Arguments:");
            System.out.println(" -a - (optional) amino acids in input FASTA (default is DNA)");
            System.out.println(" -d - (optional) print debug messages");
            System.out.println(" -m - (optional) min. number of hits in result (integer, default = )");
            System.out.println(" -M - min. sum of hit weights");
            System.out.println(" -O - order constraint");
            System.out.println(" -g - max. gap between hits to be joined");
            System.out.println(" -D - data directory with kmer-table and function-index files");
            System.out.println(" -q - ");
            System.out.println(" -");
        }
        PrintWriter pw = null;
        try {
            boolean stdout;
            if (outputFile != null) {
                pw = new PrintWriter(new FileWriter(new File(outputFile)));
                stdout = false;
            } else {
                pw = new PrintWriter(System.out);
                stdout = true;
            }
            instance.run(new File(kmerTableDir), new File(queryFastaFile), pw, stdout);
            pw.flush();
        } finally {
            if (outputFile != null) {
                pw.close();
            }
        }
    }
    
    private static void writeQueryKmer(QueryKmer item, DataOutputStream dos) throws Exception {
        dos.writeLong(item.value);
        dos.writeInt(item.hitCntId);
        dos.writeInt(item.protPos);
    }

    private static QueryKmer readQueryKmer(DataInputStream dis) throws Exception {
        QueryKmer ret = new QueryKmer();
        try {
            ret.value = dis.readLong();
            ret.hitCntId = dis.readInt();
            ret.protPos = dis.readInt();
            return ret;
        } catch (EOFException ex) {
            return null;
        }
    }

    private static void dumpQueryKmersToFile(List<QueryKmer> data, File f) throws Exception {
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(f)));
        for (QueryKmer item : data) {
            writeQueryKmer(item, dos);
        }
        dos.close();
    }
    
    private static void mergeTwoQueryKmerFiles(File f1, File f2, long numSig, File output) throws Exception {
        DataInputStream d1 = new DataInputStream(new BufferedInputStream(
                new FileInputStream(f1)));
        DataInputStream d2 = new DataInputStream(new BufferedInputStream(
                new FileInputStream(f2)));
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(output)));
        try {
            QueryKmer k1 = readQueryKmer(d1);
            QueryKmer k2 = readQueryKmer(d2);
            while (k1 != null && k2 != null) {
                long code1 = k1.value % numSig;
                long code2 = k2.value % numSig;
                int cmp = Long.compare(code1, code2);
                if (cmp == 0) {
                    cmp = Long.compare(k1.value, k2.value);
                }
                if (cmp <= 0) {
                    writeQueryKmer(k1, dos);
                    k1 = readQueryKmer(d1);
                } else {
                    writeQueryKmer(k2, dos);
                    k2 = readQueryKmer(d2);
                }
            }
            if (k1 != null) {
                writeQueryKmer(k1, dos);
            } else if (k2 != null) {
                writeQueryKmer(k2, dos);
            }
        } finally {
            d1.close();
            d2.close();
            dos.close();
        }
    }
    
    private static File mergeQueryKmerFiles(List<File> files, File tempDir, long numSig) throws Exception {
        Queue<File> filesForProcessing = new LinkedList<File>(files);
        while (filesForProcessing.size() > 1) {
            Queue<File> filesForProcessingNext = new LinkedList<File>();
            while (!filesForProcessing.isEmpty()) {
                File f1 = filesForProcessing.poll();
                if (filesForProcessing.peek() != null) {
                    File f2 = filesForProcessing.poll();
                    int chunkNum = files.size();
                    File output = new File(tempDir, "query_kmers_" + chunkNum + ".dat");
                    mergeTwoQueryKmerFiles(f1, f2, numSig, output);
                    f1.delete();
                    f2.delete();
                    files.add(output);
                    filesForProcessingNext.add(output);
                } else {
                    filesForProcessingNext.add(f1);
                }
            }
            filesForProcessing = filesForProcessingNext;
        }
        return filesForProcessing.poll();
    }
    
    public void run(File kmerTableDir, File queryFastaFile, final PrintWriter pw,
            final boolean stdout) throws Exception {
        final File tempDir = new File(".");
        File kmerTableFile = new File(kmerTableDir, "kmer.table.mem_map");
        File kmerTableGzFile = new File(kmerTableDir, kmerTableFile.getName() + ".gz");
        if (kmerTableGzFile.exists()) {
            kmerTableFile = kmerTableGzFile;
        }
        File functionIndexFile = new File(kmerTableDir, "function.index");
        File functionIndexGzFile = new File(kmerTableDir, functionIndexFile.getName() + ".gz");
        if (functionIndexGzFile.exists()) {
            functionIndexFile = functionIndexGzFile;
        }
        List<String> functionArray = loadFunctions(functionIndexFile);
        BufferedReader inputBr;
        if (queryFastaFile == null) {
            inputBr = new BufferedReader(new InputStreamReader(System.in));
        } else {
            if (queryFastaFile.getName().endsWith(".gz")) {
                inputBr = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                        new BufferedInputStream(new FileInputStream(queryFastaFile)))));
            } else {
                inputBr = new BufferedReader(new FileReader(queryFastaFile));
            }
        }
        final List<Query> queryList = new ArrayList<Query>();
        final List<HitContainer> hits = new ArrayList<HitContainer>();
        final List<QueryKmer> queryKmers = new ArrayList<QueryKmer>();
        final Map<String, Integer> queryIdToLen = new LinkedHashMap<String, Integer>();
        final KmerMemoryInfo header = new KmerMemoryInfo();
        InputStream kmerTableStream = readKmerTableHeader(kmerTableFile, header);
        final List<File> queryFiles = new ArrayList<File>();
        long t1 = System.currentTimeMillis();
        readFasta(inputBr, new FastaCallback() {
            private long loadedSoFar = 0;
            @Override
            public void nextEntry(String id, String seq, String descr) throws Exception {
                Query q = new Query();
                q.id = id;
                q.seq = seq;
                q.descr = descr;
                queryList.add(q);
                queryIdToLen.put(id, seq.length());
                loadedSoFar += seq.length();
                if (loadedSoFar > inputSizeLimit) {
                    flushChunk();
                }
            }
            @Override
            public void done() throws Exception {
                if (queryList.size() > 0) {
                    flushChunk();
                }
            }
            private void flushChunk() throws Exception {
                int chunkNum = queryFiles.size();
                prepareQueries(queryList, queryKmers, hits);
                // Let's update hash-codes of kmers and sort queries by hash-codes:
                updateHashCodeAndSort(queryKmers, header.numSigs);
                File tempFile = new File(tempDir, "query_kmers_" + chunkNum + ".dat");
                dumpQueryKmersToFile(queryKmers, tempFile);
                queryFiles.add(tempFile);
                queryKmers.clear();
                queryList.clear();
                loadedSoFar = 0;
            }
        });
        if (queryFastaFile != null) {
            inputBr.close();
        }
        File mergedQueryKmerFile = mergeQueryKmerFiles(queryFiles, tempDir, header.numSigs);
        printInfoLine("Preparation time: " + (System.currentTimeMillis() - t1) + " ms.", pw, stdout);
        long t2 = System.currentTimeMillis();
        try {
            lookup(kmerTableStream, header, mergedQueryKmerFile, hits, pw, stdout);
        } catch (Exception ex) {
            ex.printStackTrace();
            printInfoLine("Error: " + ex.getMessage(), pw, stdout);
        }
        printInfoLine("Lookup time: " + (System.currentTimeMillis() - t2) + " ms.", pw, stdout);
        long t3 = System.currentTimeMillis();
        Map<HitContainerKey, HitContainer> hitCnts = 
                new LinkedHashMap<HitContainerKey, HitContainer>();
        for (HitContainer cnt : hits) {
            hitCnts.put(cnt.key, cnt);
        }
        for (String queryId : queryIdToLen.keySet()) {
            int seqLen = queryIdToLen.get(queryId);
            if (aa) {
                processAASeq(queryId, seqLen, hitCnts, functionArray, pw);
            } else {
                processSeq(queryId, seqLen, hitCnts, functionArray, pw);
            }
            pw.flush();
        }
        printInfoLine("Grouping time: " + (System.currentTimeMillis() - t3) + " ms.", pw, stdout);
    }
    
    private void printInfoLine(String message, PrintWriter pw, boolean stdout) {
        if (debug) {
            pw.println(message);
        }
        if (!stdout) {
            System.out.println(message);
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
        hitCnt.id = hitCnts.size();
        hitCnts.add(hitCnt);
        for (int i = 0; i < pIseq.length - K; i++) {
            long value = encodedKmer(pIseq, i);
            if (value < 0)
                continue;
            QueryKmer qk = new QueryKmer();
            qk.value = value;
            qk.protPos = i;
            qk.hitCntId = hitCnt.id;
            queryList.add(qk);
        }
    }
    
    private InputStream readKmerTableHeader(File kmerTableFile, 
            KmerMemoryInfo ret) throws Exception {
        InputStream is;
        if (kmerTableFile.getName().endsWith(".gz")) {
            is = new BufferedInputStream(new GZIPInputStream(new BufferedInputStream(
                    new FileInputStream(kmerTableFile))));
        } else {
            is = new BufferedInputStream(new FileInputStream(kmerTableFile));
        }
        long numSigs = readLongLE(is);
        long entrySize = readLongLE(is);
        long version = readLongLE(is);
        if (ret != null) {
            ret.numSigs = numSigs;
            ret.entrySize = entrySize;
            ret.version = version;
        }
        return is;
    }
    
    private void lookup(InputStream kmerTableStream, KmerMemoryInfo header, 
            File sortedQueryKmerFile, List<HitContainer> hitCnts, 
            PrintWriter pw, boolean stdout) throws Exception {
        InputStream is = kmerTableStream;
        long numSigs = header.numSigs;
        long entrySize = header.entrySize;
        long version = header.version;
        if (debug) {
            pw.println("Kmer-table info: numSigs=" + numSigs + ", entrySize=" + entrySize + 
                    ", version=" + version);
        }
        
        long t1 = System.currentTimeMillis();
        int kmersFound = 0;
        int posCount = 0;
        DataInputStream queryStream = new DataInputStream(new BufferedInputStream(
                new FileInputStream(sortedQueryKmerFile)));
        try {
            // Now we go along hash table and along query hash-codes
            long curHashCode = 0;
            QueryKmer curKmer = readQueryKmer(queryStream);
            Map<Long, List<QueryKmer>> inProgress = new HashMap<Long, List<QueryKmer>>();
            int fraction = 0;
            while (curKmer != null || inProgress.size() > 0) {
                long neededHashCode = curHashCode;
                if (inProgress.size() == 0) {
                    // Update next hash-code
                    QueryKmer qk = curKmer;
                    neededHashCode = qk.value % numSigs;
                    List<QueryKmer> list = new ArrayList<QueryKmer>(5);
                    list.add(qk);
                    inProgress.put(qk.value, list);
                    curKmer = readQueryKmer(queryStream);
                }
                // Let's push all queries with necessary hash-code into progress state
                while (curKmer != null) {
                    QueryKmer qk = curKmer;
                    if ((qk.value % numSigs) != neededHashCode) {
                        break;
                    }
                    if (inProgress.containsKey(qk.value)) {
                        inProgress.get(qk.value).add(qk);
                    } else {
                        List<QueryKmer> list = new ArrayList<QueryKmer>(3);
                        list.add(qk);
                        inProgress.put(qk.value, list);
                    }
                    curKmer = readQueryKmer(queryStream);
                }
                // Let's position kmer-table stream to necessary hash-code
                if (neededHashCode > curHashCode) {
                    skipBytesFully(is, entrySize * (long)(neededHashCode - curHashCode));
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
                        for (QueryKmer qk : inProgress.remove(whichKmer)) {
                            Hit hit = new Hit();
                            hit.from0InProt = qk.protPos;
                            hit.oI = otuIndex;
                            hit.avgOffFromEnd = avgFromEnd;
                            hit.fI = functionIndex;
                            hit.functionWt = functionWt;
                            hitCnts.get(qk.hitCntId).hits.add(hit);
                            posCount++;
                        }
                    }
                }
                curHashCode++;
                int newFraction = (int)(100.0 * ((double)curHashCode / (double)numSigs));
                if (newFraction != fraction) {
                    fraction = newFraction;
                    printInfoLine("Processed: " + (fraction) + "%, time=" +
                            (System.currentTimeMillis() - t1) + " ms., found-so-far=" + 
                            kmersFound, pw, stdout);
                }
            }
        } finally {
            is.close();
            queryStream.close();
        }
        if (debug) {
            pw.println("Kmers found: " + kmersFound + " (pos-count=" + posCount + ")");
        }
    }

    private static void skipBytesFully(InputStream is, long bytesToSkip) throws IOException {
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
    }

    private void prepareQueries(List<Query> input, List<QueryKmer> queryKmers,
            List<HitContainer> hitCnts) {
        for (Query q : input) {
            String id = q.id;
            char[] seq = q.seq.toCharArray();
            if (aa) {
                byte[] pIseq = new byte[seq.length];
                for (int i = 0; i < seq.length; i++) {
                    pIseq[i] = toAminoAcidOff(seq[i]);
                }
                addKmers(id, '+', 0, seq, pIseq, queryKmers, hitCnts);
            } else {
                int len = seq.length / 3 + 1;
                char[] pseq = new char[len];
                byte[] pIseq = new byte[len];
                for (int frame = 0; frame < 3; frame++) {
                    translate(seq, frame, pseq, pIseq);
                    addKmers(id, '+', frame, pseq, pIseq, queryKmers, hitCnts);
                }
                char[] complSeq = revComp(seq);
                for (int frame = 0; frame < 3; frame++) {
                    translate(complSeq, frame, pseq, pIseq);
                    addKmers(id, '-', frame, pseq, pIseq, queryKmers, hitCnts);
                }
            }
        }
    }

    private static void updateHashCodeAndSort(List<QueryKmer> values, final long numSigs) {
        // Sort query kmers by hash-codes (to be able to process them in parrallel with
        // processing kmer-table file (sorted joining).
        Collections.sort(values, new Comparator<QueryKmer>() {
            @Override
            public int compare(QueryKmer o1, QueryKmer o2) {
                long hashCode1 = o1.value % numSigs;
                long hashCode2 = o2.value % numSigs;
                int ret = Long.compare(hashCode1, hashCode2);
                if (ret == 0) {
                    ret = Long.compare(o1.value, o2.value);
                }
                return ret;
            }
        });
    }

    private static int readIntLE(InputStream is) throws IOException {
        int ch1 = is.read();
        int ch2 = is.read();
        int ch3 = is.read();
        int ch4 = is.read();
        if ((ch1 | ch2 | ch3 | ch4) < 0)
            throw new EOFException();
        return ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
    }

    private static long readLongLE(InputStream is) throws IOException {
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

    private static float readFloatLE(InputStream is) throws IOException {
        return Float.intBitsToFloat(readIntLE(is));
    }
    
    private static void readFasta(BufferedReader br, FastaCallback ret) throws Exception {
        String str1 = null;
        while (true) {
            String protName = null;
            String protSeq = null;
            String protDescr = null;
            try {
                if(str1 == null)
                    str1 = br.readLine();
                for(;;) {
                    if(str1 == null)
                        break;
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
                if (protName == null) {
                    break;
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
                ret.nextEntry(protName, protSeq, protDescr);
            } catch (IllegalStateException ex) {
                throw ex;
            } catch (Exception ex) {
                throw new IllegalStateException(ex);
            }
        }
        ret.done();
    }
    
    public static class KmerMemoryInfo {
        public long numSigs;
        public long entrySize;
        public long version;
    }
    
    public static class Query {
        public String id;
        public String descr;
        public String seq;
    }
    
    public static class QueryKmer {
        public long value;
        public int hitCntId;
        public int protPos;
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
        public int avgOffFromEnd;    /* average offset from the end */
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
        public int id;
        public List<Hit> hits;
    }
    
    public static interface FastaCallback {
        public void nextEntry(String id, String seq, String descr) throws Exception;
        public void done() throws Exception;
    }
    
    /*public static interface QueryKmerSorter {
        public void addKmer(QueryKmer qk) throws Exception;
        public void finalizeSorting() throws Exception;
        public QueryKmer loadNext() throws Exception;
    }*/
}