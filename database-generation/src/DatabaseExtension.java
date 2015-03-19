import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;

public class DatabaseExtension {	
	private static Charset UTF8 = Charset.forName("UTF-8");
	private String[][][] originalLines;
	private String[][] destinationLines;
		
	private FileWriter[] writers;
	private FileWriter writer;
	private int variants;
	private int currentWriter = -1;
	private int lastLine = 0;
	
	private long startingTime = 0;
	public long deltaTime = 0;
	
	public DatabaseExtension(){
		originalLines = new String[Main.THREADS][Main.LINES_BY_THREAD][8];
		for(int i=0; i < Main.THREADS; i++)
			for(int j=0; j < Main.LINES_BY_THREAD; j++)
				for(int k=0; k < 8; k++)
					originalLines[i][j][k] = "";
		
		destinationLines = new String[Main.THREADS][Main.LINES_BY_THREAD];	
		for(int i=0; i < Main.THREADS; i++)
			for(int j=0; j < Main.LINES_BY_THREAD; j++)
				destinationLines[i][j] = "";		
	}
			

	/**
	 * Read the vcf file, then call FilePart to analyze & modify the data, and save them
	 * @return
	 */
	public boolean doJob() {
        String line;
        String[] elements = new String[1];
        int i=0, j=0, destinationThread=-1;
        
		// We open the destination file (creating/truncating)
		PrintWriter writ;
		try {
			if(Main.ONE_BIG_FILE == true) {
				writ = new PrintWriter(Main.destination, "UTF-8");
				writ.close();
			} else {
				for(i=0; i < Main.MAX_FILES; i++) {
					writ = new PrintWriter(Main.destination.concat(".").concat(String.valueOf(i)), "UTF-8");
					writ.close();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//We keep a pointer to the destination file(s)
		if(Main.ONE_BIG_FILE == true){
			try {
				writer = new FileWriter(Main.destination, true);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else {
			writers = new FileWriter[Main.MAX_FILES];
			for( i=0; i < Main.MAX_FILES; i++) {
				try {
					writers[i] = new FileWriter(Main.destination.concat(".".concat(String.valueOf(i))), true);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
				
		//We check if the origin file exists at least
		File f = new File(Main.origin);
		if(!f.exists())
			return false;

		startingTime = System.nanoTime();
		//We read the file line by line
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(Main.origin), UTF8));  
			i=0;
		    while ((line = br.readLine()) != null) {
		    	
		    	//We check which thread we will use
		    	destinationThread = (int) Math.floor(i/Main.LINES_BY_THREAD);
		    	
		    	//We check if we have the header or some interesting information
		    	elements = line.trim().split("\t");
    			this.lastLine++;
		    	
		    	if(elements.length >= 8) {		    		
		    		if(!"INFO".equals(elements[7])) {
		    			j = i - destinationThread*Main.LINES_BY_THREAD;
		    			originalLines[destinationThread][j] = elements;
		    			i++;		    			
		    		}
		    		else {
		    			this.generateFirstLine(line);
		    		}		    		
		    	} else {
		    		this.writeHeader(line);
		    	}
		    	
		    	//If we have enough data, we prepare the threads and reinitialize the counter
		    	if(i >= Main.LINES_BY_THREAD * Main.THREADS) {
		    		this.launchThreads(destinationThread+1);
		    		i = 0;
		    	}
		    }
		    
		    if(i > 0)
		    	this.launchThreads(destinationThread+1);		    
		    
		    br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return true;
	}
	
	/**
	 * Prepare and launch the threads
	 * @param lastThread the last thread with some information
	 */
	public void launchThreads(int lastThread) {
		ArrayList<Thread> fileparts = new ArrayList<Thread>();

		for(currentWriter = 0; currentWriter < Main.MAX_FILES; currentWriter++) {
			if(Main.ONE_BIG_FILE == true)
				currentWriter = -1;

			for(int i=0; i < lastThread; i++){
				FilePart filepart = new FilePart(i, this, this.originalLines[i]);
				
				Thread t = new Thread(filepart);
				fileparts.add(t);
			    t.start();
			}
			
	    	//We wait for all the threads to finish their processing
			try {
				for(Thread t : fileparts)
					t.join();
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
					
	    	//We save the data
			try {
				writeVariantsToDisk();
			} catch (IOException e) {
				e.printStackTrace();
			}		
			
			if(Main.ONE_BIG_FILE == true)
				break;
		}		
		this.displayExecutionTime();
		
		//We reinitialise the data
		originalLines = new String[Main.THREADS][Main.LINES_BY_THREAD][8];
		for(int i=0; i < Main.THREADS; i++)
			for(int j=0; j < Main.LINES_BY_THREAD; j++)
				for(int k=0; k < 8; k++)
					originalLines[i][j][k] = "";
	}
	
	/**
	 * Write the variants to disk
	 * @throws IOException 
	 */
	public void writeVariantsToDisk() throws IOException {
		for(int thread=0; thread < this.destinationLines.length; thread++) {	
			for(int j=0; j < this.destinationLines[thread].length; j++) {

				if(!this.destinationLines[thread][j].isEmpty() && Main.THREADS > 1){//We can have no variant on some lines
					if(currentWriter == -1)
						writer.write(destinationLines[thread][j]);
					else
						writers[currentWriter].write(destinationLines[thread][j]);
					variants++;
				}
				this.destinationLines[thread][j] = "";
			}
		}
	}
	
	/**
	 * Saves a specific line in memory from a thread
	 */
	public void setDestinationLine(int thread, int lineId, String line) {
		if(Main.THREADS == 1) {
			//If we have only one thread, we can directly write into the file
			try {
				if(Main.ONE_BIG_FILE == true)
					writer.write(line);
				else 
					writers[currentWriter].write(line);
				variants++;
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else
			this.destinationLines[thread][lineId] = line; 
	}
	
	/**
	 * Generate the first line of interesting information and write it to the file
	 * @param line the initial line where we need to add the user ids
	 */
	private void generateFirstLine(String line) {
		
		//String.format() is very slow but it's just for 1 line.
		for(int i=1; i <= Main.SAMPLES; i++)
			line += "\tNA"+String.format("%05d", i);
		
		try {
			if(Main.ONE_BIG_FILE == true)
				writer.write(line.concat("\r\n"));
			else {
				for(int i=0; i < Main.MAX_FILES; i++)
					writers[i].write(line.concat("\r\n"));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * We manage the header directly here as FilePart is optimized for variant lines.
	 */
	private void writeHeader(String line) {
		try {
			if(Main.ONE_BIG_FILE == true)
				writer.write(line.concat("\r\n"));
			else {
				for(int i=0; i < Main.MAX_FILES; i++)
					writers[i].write(line.concat("\r\n"));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public int getVariants() {
		return variants;
	}
	
	public void displayExecutionTime() {
		if(deltaTime > 0)
			System.out.println("Total execution time: "+String.valueOf(Math.round((System.nanoTime() - startingTime)/1000000000))+"s with subtotal time: "+String.valueOf(Math.round(deltaTime/1000000000))+"s. Progression: "+String.valueOf(Math.round(this.lastLine*10000.0/9462924)/100.0)+"%.");
		else
			System.out.println("Total execution time: "+String.valueOf(Math.round((System.nanoTime() - startingTime)/1000000000))+"s. Progression: "+String.valueOf(Math.round(this.lastLine*10000.0/9462924)/100.0)+"%.");
		
	}
}
