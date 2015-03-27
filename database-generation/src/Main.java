import java.io.File;

import org.apache.commons.lang3.StringUtils;


public class Main {
	public static int THREADS = 5;
	public static int SAMPLES = 1000;
	public static boolean ONE_BIG_FILE = true;//Indicates if we create only one big file at the output, or 1 file/sample
	
	public static String origin = "../../Public_Database/Public_Database_v3";
	public static String destination = "../../Public_Database/Public_Database_v3_"+Main.SAMPLES;

	//Computed automatically, no need to modify it
	public static int MAX_FILES = 250;
	public static int LINES_BY_THREAD = 10;
	public static int TOTAL_SAMPLES = 60706;
	
	public static long rand = 0;
	
	public static void main(String[] args) {

		//Main.stringSplit();
		//Main.stringConcatenation();
		//System.exit(0);
		
		//We analyze the received parameters
		if(args.length != 5*2) {
			System.err.println("Invalid number of parameters. You have to give --o [path-to-database] --d [path-to-future-vcf] --s [number-of-samples] --f [boolean-big-file-or-small-files] --t [number-of-threads]");
			System.exit(0);
		} else {
			boolean error = false;
			int differentArguments = 0;
			
			//We check each parameter
			try {
				for(int i=0; i < args.length; i += 2) {
					if("--o".equals(args[i])) {
						
						if(!new File(args[i+1]).exists()) {
							System.err.println("Invalid path to database given. File not found");
							error = true;
						} else {
							origin = args[i+1];
						}		
						differentArguments++;
						
					} else if("--d".equals(args[i])) {
						
						destination = args[i+1];
						differentArguments++;
						
					} else if("--s".equals(args[i])) {
						try {
							Main.SAMPLES = Integer.valueOf(args[i+1]);
						} catch (Exception e) {
							System.err.println("Invalid number of samples to generate.");
							error = true;
						}
						
						if(Main.SAMPLES <= 0) {
							System.err.println("Invalid number of samples to generate.");	
							error = true;
						} else if(Main.SAMPLES > 10000) {
							System.out.println("You have asked to generate a vcf for "+String.valueOf(Main.SAMPLES)+" samples. I hope for you that you have enough space to save it...");
						}
						
						differentArguments++;
					} else if("--f".equals(args[i])) {
						if("1".equals(args[i+1]) || "true".equals(args[i+1].toLowerCase())) {
							Main.ONE_BIG_FILE = true;
						} else if("0".equals(args[i+1]) || "false".equals(args[i+1].toLowerCase())) {
							Main.ONE_BIG_FILE = false;
						} else {
							System.out.println("The parameter --f is asking for a boolean so, you have to give (1, true, 0 or false)");
							error = true;
						}
						
						differentArguments++;
					} else if("--t".equals(args[i])) {
						try {
							Main.THREADS = Integer.valueOf(args[i+1]);
						} catch (Exception e) {
							System.err.println("Invalid number of threads to use.");
							error = true;
						}
						
						if(Main.THREADS <= 0) {
							System.err.println("Invalid number of threads to use.");	
							error = true;
						} else if(Main.THREADS >= 100) {
							System.out.println("You have asked to use "+String.valueOf(Main.THREADS)+" threads... Are you crazy or do you come from the future? Nevermind, do as you please...");
						}
						
						differentArguments++;
					}
				}
			} catch (Exception e) {
				System.err.println("Invalid parameters. You have to give --o [path-to-database] --d [path-to-future-vcf] --s [number-of-samples] --f [boolean-big-file-or-small-files] --t [number-of-threads]");
				System.out.println(e.toString());
				System.exit(0);
			}
			
			if(error || differentArguments != 5) {
				System.err.println("Invalid parameters. You have to give --o [path-to-database] --d [path-to-future-vcf] --s [number-of-samples] --f [boolean-big-file-or-small-files] --t [number-of-threads]");
				System.exit(0);
			}
		}

		//We compute the number of lines each thread can use, based on the average length size
		int averageLength = 2*1024 + Main.SAMPLES * 14;
		Main.LINES_BY_THREAD = (int)(0.1*(Runtime.getRuntime().freeMemory()/(Main.THREADS * averageLength)));		
		destination += "_s"+String.valueOf(Main.SAMPLES);
		
		if(ONE_BIG_FILE == false) {
			MAX_FILES = SAMPLES;
			SAMPLES = 1;
		} else {
			destination += ".vcf";
		}

		//We launch the job
		System.out.println("We launch the job, please be patient...");
		DatabaseExtension databaseExtension = new DatabaseExtension();
		boolean res = databaseExtension.doJob();
		
		if(res)
			System.out.println("Job done with "+String.valueOf(databaseExtension.getVariants())+" variants.");
		else
			System.out.println("Job not done.");
		
		databaseExtension.displayExecutionTime();
	}
		
	public static void stringConcatenation() {
		String[] tab = new String[20000000];
		String[] tb = new String[tab.length];
		String text;
		int i, j, lines=900, samples=60000;
		long start;
		
		for(i=0; i < tab.length; i++)
			tab[i] = String.valueOf((int) (Math.random()*i));

		start = System.currentTimeMillis();
		for(i=0; i < 100000000; i++)
			Math.random();
		System.out.println("Method 0a (random): Total time for 100 000 000 iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		start = System.currentTimeMillis();
		Main.XORShift64Init(System.nanoTime());
		for(i=0; i < 100000000; i++)
			Main.random01();// Interesting to notice: xorshift is 10x faster than random() on my i7 2.8GHz, but slower than random() on a i5 3.5GHz... 
		System.out.println("Method 0b (xorshift random): Total time for 100 000 000 iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		start = System.currentTimeMillis();
		for(i=0; i < tab.length; i++)
			tb[i] = tab[i];
		System.out.println("Method 0c (assignation, not concatenation): Total time for "+tab.length+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		start = System.currentTimeMillis();
		for(j=0; j < 9; j++) {
			text = "";
			for(i=0; i < samples; i++)
				text += tab[i];
		}
		System.out.println("Method 1: Total time for "+(9*samples)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		start = System.currentTimeMillis();
		for(j=0; j < 9; j++) {
			text = "";
			for(i=0; i < samples; i++)
				text = text.concat(tab[i]);
		}
		System.out.println("Method 2: Total time for "+(9*samples)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		
		start = System.currentTimeMillis();
		StringBuffer t;
		for(j=0; j < lines; j++) {
			t = new StringBuffer();
			for(i=0; i < samples; i++)
				t = t.append(tab[i]);
			text = t.toString();
		}
		System.out.println("Method 3: Total time for "+(lines*samples)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");
		
		// Something weird... If I use this below, it is slower than method 3. But in the real code, it's almost the same, even the opposite when we have 20.000 samples...
		start = System.currentTimeMillis();
		for(j=0; j < 9; j++) {
			for(i=1; i < samples; i++)
				tb[i] = new StringBuffer(tab[i]).toString();
			text = StringUtils.join(tb);
		}
		text = StringUtils.join(tb);
		System.out.println("Method 4: Total time for "+(9*samples)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");		
	}

	@SuppressWarnings("unused")
	public static void stringSplit() {
		String line = "1	13528	.	C	G,T	1771.54	VQSRTrancheSNP99.60to99.80	AC=22,11;AC_AFR=13,0;AC_AMR=1,0;AC_Adj=14,9;AC_EAS=0,0;AC_FIN=0,0;AC_Het=14,9,0;AC_Hom=0,0;AC_NFE=0,2;AC_OTH=0,0;AC_SAS=0,7;AF=6.158e-04,3.079e-04;AN=35726;AN_AFR=492;AN_AMR=124;AN_Adj=11114;AN_EAS=154;AN_FIN=8;AN_NFE=3122;AN_OTH=126;AN_SAS=7088;BaseQRankSum=1.23;ClippingRankSum=0.056;DP=152693;FS=0.000;GQ_MEAN=14.54;GQ_STDDEV=16.53;Het_AFR=13,0,0;Het_AMR=1,0,0;Het_EAS=0,0,0;Het_FIN=0,0,0;Het_NFE=0,2,0;Het_OTH=0,0,0;Het_SAS=0,7,0;Hom_AFR=0,0;Hom_AMR=0,0;Hom_EAS=0,0;Hom_FIN=0,0;Hom_NFE=0,0;Hom_OTH=0,0;Hom_SAS=0,0;InbreedingCoeff=0.0557;MQ=31.08;MQ0=0;MQRankSum=-5.410e-01;NCC=67387;QD=1.91;ReadPosRankSum=0.206;VQSLOD=-2.705e+00;culprit=MQ;DP_HIST=10649|1543|748|1347|2650|650|180|54|18|12|8|3|0|0|1|0|0|0|0|0,2|6|2|1|4|0|4|1|0|0|2|0|0|0|0|0|0|0|0|0,1|0|0|0|1|1|3|0|1|1|1|0|0|0|1|0|0|0|0|0;GQ_HIST=345|11296|88|59|3313|547|382|61|12|4|5|7|1501|198|18|16|1|0|1|9,0|0|1|0|1|0|3|1|0|1|2|0|1|2|0|1|1|0|1|7,0|1|0|0|1|1|0|0|1|0|0|1|1|1|1|0|0|0|0|2;CSQ=T|ENSG00000223972|ENST00000456328|Transcript|non_coding_exon_variant&nc_transcript_variant|776||||||||2|3/3|||||||1||YES|DDX11L1|HGNC||||processed_transcript||||ENST00000456328.2:n.776C>T||||||||||,G|ENSG00000223972|ENST00000456328|Transcript|non_coding_exon_variant&nc_transcript_variant|776||||||||1|3/3|||||||1||YES|DDX11L1|HGNC||||processed_transcript||||ENST00000456328.2:n.776C>G||||||||||,T|ENSG00000223972|ENST00000450305|Transcript|non_coding_exon_variant&nc_transcript_variant|490||||||||2|6/6|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000450305.2:n.490C>T||||||||||,G|ENSG00000223972|ENST00000450305|Transcript|non_coding_exon_variant&nc_transcript_variant|490||||||||1|6/6|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000450305.2:n.490C>G||||||||||,T|ENSG00000223972|ENST00000515242|Transcript|non_coding_exon_variant&nc_transcript_variant|769||||||||2|3/3|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000515242.2:n.769C>T||||||||||,G|ENSG00000223972|ENST00000515242|Transcript|non_coding_exon_variant&nc_transcript_variant|769||||||||1|3/3|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000515242.2:n.769C>G||||||||||,T|ENSG00000223972|ENST00000518655|Transcript|non_coding_exon_variant&nc_transcript_variant|607||||||||2|3/4|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000518655.2:n.607C>T||||||||||,G|ENSG00000223972|ENST00000518655|Transcript|non_coding_exon_variant&nc_transcript_variant|607||||||||1|3/4|||||||1|||DDX11L1|HGNC||||transcribed_unprocessed_pseudogene||||ENST00000518655.2:n.607C>G||||||||||,T||ENSR00000528767|RegulatoryFeature|regulatory_region_variant|||||||||2||||||||||||||||||||||||||||||,G||ENSR00000528767|RegulatoryFeature|regulatory_region_variant|||||||||1||||||||||||||||||||||||||||||";
		String tmp[], afValue;
		int lines = 9, samples=60000, i, j, startIndex, endIndex;
		long start;
		
		start = System.currentTimeMillis();
		for(i=0; i < lines; i++)
			for(j=0; j < samples; j++)
				tmp = line.split("	");
		System.out.println("Method 1a: Total time for "+String.valueOf(lines)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");	
		
		start = System.currentTimeMillis();
		for(i=0; i < lines; i++) {
			for(j=0; j < samples; j++) {
				tmp = line.split("AF=");
				tmp = tmp[1].split(";");
				afValue = tmp[0];
			}
		}
		System.out.println("Method 2a: Total time for "+String.valueOf(lines)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");	
			
		start = System.currentTimeMillis();
		for(i=0; i < lines; i++) {
			for(j=0; j < samples; j++) {
				startIndex = line.indexOf("AF=") + 3;
				endIndex = line.indexOf(";", startIndex);
				afValue = line.substring(startIndex, endIndex);
			}
		}
		System.out.println("Method 2b: Total time for "+String.valueOf(lines)+" iterations: "+((System.currentTimeMillis()-start)/1000.0)+"s.");	
	}
	
	/**
	 * Initialize a random generator better and faster (> 10x) than Math.random()
	 * @see http://demesos.blogspot.be/2011/09/replacing-java-random-generator.html
	 */
	public static void XORShift64Init(long seed) {
		Main.rand = seed==0 ? 0xdeadbeef : seed;
	}

	/**
	 * Return a random value better and faster (> 10x) than Math.random()
	 * @see http://demesos.blogspot.be/2011/09/replacing-java-random-generator.html
	 * @see http://en.wikipedia.org/wiki/Xorshift
	 */
	public static double random01() {
		Main.rand ^= (Main.rand << 21); 
		Main.rand ^= (Main.rand >>> 35);
		Main.rand ^= (Main.rand << 4);
		
		return (Main.rand*1.0 / Long.MAX_VALUE)/2.0 + 0.5;
	}
}
