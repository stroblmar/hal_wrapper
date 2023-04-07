package onLatticeCA_jar;

import HAL.lib.CommandLine;
import java.io.File;
import static HAL.Util.LoadState;

// ------------------------------------------------------------------------------------------------------------
// Main "Wrapper" class that manages input/output and the simulation.
// ------------------------------------------------------------------------------------------------------------
@CommandLine.Command(name = "Game Assay Simulator",
        mixinStandardHelpOptions = true,
        showDefaultValues = true,
        description = "A 2D on-lattice agent-based model that simulates 2 cell populations growing together in a domain under a given drug schedule.")
public class OnLatticeModel implements Runnable {
    /* ========================================================================
     * --- Parameters ----
     * ======================================================================== */
    // Parse variables from command line using piclocli (http://picocli.info)
    OnLatticeGrid model = new OnLatticeGrid(); // Temporary placeholder so can load defaults from the base class
    // ------------------------- Config -------------------------
    @CommandLine.Option(names = { "-h", "--help" }, usageHelp = true, description = "Display a help message and exit.") 
    boolean helpRequested;
    @CommandLine.Option(names = { "-s", "--seed"}, description="Random number seed.") 
    int seed = model.seed;
    @CommandLine.Option(names = { "--headless" }, description = "Whether or not to run in headless mode (used for running on cluster).") 
    boolean headless = true;

    // ------------------------- Experimental Setup -------------------------
    @CommandLine.Option(names = { "-xDim", "--xDim"}, description="x-dimension of domain (in lattice sites)") 
    int xDim = model.xDim;
    @CommandLine.Option(names = { "-yDim", "--yDim"}, description="y-dimension of domain (in lattice sites)") 
    int yDim = model.yDim;
    @CommandLine.Option(names = { "-n", "--nReplicates"}, description="Number of replicates.") 
    int nReplicates = model.nReplicates;
    @CommandLine.Option(names = { "-N0", "--initialSize"}, description="Total initial cell number (absolute)") 
    double initialSize = model.initialSize;
    @CommandLine.Option(names = { "-n0", "--initialSizeProp"}, description="Initial cell density relative to (physical) carrying capacity") 
    double initialSizeProp = model.initialSizeProp;
    @CommandLine.Option(names = { "-fR", "--rFrac"}, description="Initial resistance fraction in [0,1]") 
    double rFrac = model.rFrac;    
    @CommandLine.Option(names = { "-tEnd", "--tEnd"}, description="End time in days.") 
    double tEnd = model.tEnd;
    @CommandLine.Option(names = { "-dt", "--dt"}, description="Time step in days.") 
    double dt = model.dt;
    @CommandLine.Option(names = { "-treatmentScheduleList", "--treatmentScheduleList"}, description="Treatment schedule in format {{tStart, tEnd, drugConcentration}}.") 
    String treatmentScheduleList_string;

    // ------------------------- Cell Properties -------------------------
    @CommandLine.Option(names = { "-r_S", "--divisionRate_S"}, description="Proliferation rate of sensitive cells in d^-1. Proliferation will be attempted at this rate.") 
    double divisionRate_S = model.divisionRate_S;
    @CommandLine.Option(names = { "-r_R", "--divisionRate_R"}, description="Proliferation rate of sensitive cells in d^-1. Proliferation will be attempted at this rate.") 
    double divisionRate_R = model.divisionRate_R;
    @CommandLine.Option(names = { "-m_S", "--movementRate_S"}, description="Movement rate of sensitive cells in d^-1. Movement will be attempted at this rate.") 
    double movementRate_S = model.movementRate_S;
    @CommandLine.Option(names = { "-m_R", "--movementRate_R"}, description="Movement rate of resistant cells in d^-1. Movement will be attempted at this rate.") 
    double movementRate_R = model.movementRate_R;
    @CommandLine.Option(names = { "-d_S", "--deathRate_S"}, description="Natural death rate of sensitive cells in d^-1. Mean life span of a cell will be 1/k.") 
    double deathRate_S = model.deathRate_S;
    @CommandLine.Option(names = { "-d_R", "--deathRate_R"}, description="Natural death rate of resistant cells in d^-1. Mean life span of a cell will be 1/k.") 
    double deathRate_R = model.deathRate_R;
    @CommandLine.Option(names = { "-d", "--drugKillProportion"}, description="Drug induced death rate in d^-1.") 
    double drugKillProportion = model.drugKillProportion;

    // ------------------------- Output - Text -------------------------
    @CommandLine.Option(names = { "-v", "--verboseLevel"}, description="Verbosity (0-2).") 
    int verboseLevel = model.verboseLevel;
    @CommandLine.Option(names = { "-printfreq", "--printFrequency"}, description="Frequency at which output is printed to the screen.") 
    double printFrequency = model.printFrequency;
    @CommandLine.Option(names = { "-logf", "--logCellCountFrequency"}, description="Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging.") 
    double logCellCountFrequency = model.logCellCountFrequency;
    // @CommandLine.Option(names = { "--extraSimulationInfo"}, description="Additional columns to be printed into the log file. Useful when want to annotate experiment to facilitate later analysis.") 
    // double[] extraSimulationInfo = model.extraSimulationInfo;
    // @CommandLine.Option(names = { "--extraSimulationInfoNames"}, description="Column headers to be used for additional columns to be printed into the log file.") 
    // String[] extraSimulationInfoNames = model.extraSimulationInfoNames;
    @CommandLine.Option(names = { "--outDir"}, description="Directory which to save output files to.") 
    String outDir = "./tmp/";
    
    // ------------------------- Output - Visualisation -------------------------
    @CommandLine.Option(names = { "-visualiseB", "--visualiseB"}, description="Whether or not to show visualization.")
    Boolean visualiseB = model.visualiseB;
    @CommandLine.Option(names = { "--scaleFactor"}, description="Scale factor used to display grid.") 
    int scaleFactor = model.scaleFactor;
    @CommandLine.Option(names = { "--pause"}, description="Pause betwween time steps to slow simulation.") 
    int pause = model.pause;
    @CommandLine.Option(names = { "--imageFrequency"}, description="Frequency at which an image of the tumour is saved. Negative number turns it off.") 
    int imageFrequency = model.imageFrequency;
    @CommandLine.Option(names = { "--imageOutDir"}, description="Directory which to save images to.") 
    String imageOutDir = model.imageOutDir;

    // ------------------------- Output - Model File -------------------------
    @CommandLine.Option(names = { "--saveModelState"}, description="Whether or not to save the model object at the end of the simulation.") 
    Boolean saveModelState = model.saveModelState;
    // @CommandLine.Option(names = { "--fromScratch"}, description="Whether or not to run simulation from scratch or from a model file.") 
    // Boolean fromScratch = model.fromScratch; // #TODO: Can I clean this up?
    @CommandLine.Option(names = { "--savedModelFileName"}, description="Name of model file to load when continuing a previous run.") 
    String savedModelFileName = model.savedModelFileName;

    // ------------------------------------------------------------------------------------------------------------
    /*
     * Main runner function
     */
    public void run(){
        // ------------------------- 1. Setup -------------------------
        assert !helpRequested; // Print help message if requested
        if (headless) {System.setProperty("java.awt.headless", "true");} // Make headless if necessary
        OnLatticeGrid myModel;
        double[] paramArr;
        String currSavedModelFileName;
        String outFName;
        Boolean fromScratch = savedModelFileName==null;
        new File(outDir).mkdirs(); // Set up environment

        // Parse treatment schedule. If none set, simulate growth until the end of the simulation.
        double[][] treatmentScheduleList;
        if (treatmentScheduleList_string == null) {
            treatmentScheduleList = new double[][] {{0, tEnd, 0}};
        } else {
            // Parse command line input
            String[] tmpList = treatmentScheduleList_string.split("\\[");
            String[] processed_justNumbers;
            treatmentScheduleList = new double[tmpList.length-2][3];
            for (int k=2; k<tmpList.length; k++) {
                processed_justNumbers = tmpList[k].split("\\]");
                processed_justNumbers = processed_justNumbers[0].split(",");
                for (int kk=0; kk<processed_justNumbers.length; kk++) {
                    treatmentScheduleList[k-2][kk] = Double.parseDouble(processed_justNumbers[kk]);
                }
            }
        }

        // Parse initial size information. If fraction provided use that. Otherwise use absolute numbers
        if (initialSizeProp > 0) {
            initialSize = (int) (initialSizeProp * xDim * yDim);
        }
        int[] initPopSizeArr = new int[]{(int) Math.floor(initialSize * (1 - rFrac)), (int) Math.ceil(initialSize * rFrac)};

        // Prepare storage and indexing for each replicate
        int[] replicateIdList = new int[] {seed};
        int replicateId;
        if (nReplicates>1) {
            replicateIdList = new int[nReplicates];
            for (int i=0; i<nReplicates; i++) {replicateIdList[i]=i;}
        }

        // ------------------------------------------------------------------------
        // Main simulation loop.
        // ------------------------------------------------------------------------
        for (int replicateIdx = 0; replicateIdx < nReplicates; replicateIdx++) {
            replicateId = replicateIdList[replicateIdx];
            currSavedModelFileName = savedModelFileName==null? outDir + "RepId_" + replicateId + ".bin": savedModelFileName;

            // Models can be created either from scratch, or be loaded from a previous run
            if (fromScratch) {
                // Set up the simulation
                paramArr = new double[]{divisionRate_S, divisionRate_R, movementRate_S, movementRate_R,
                        deathRate_S, deathRate_R, drugKillProportion};
                myModel = new OnLatticeGrid(xDim, yDim, paramArr, dt);
                myModel.tEnd = tEnd;
                myModel.visualiseB = visualiseB;
                myModel.verboseLevel = verboseLevel;
                myModel.pause = pause;

                // Set the random number seed. Behaviour depends no whether this is a single run or part of a series of nReplicate runs. By default will assign every replicate the value ```seed=replicateId```
                if (seed != -1) {
                    if (nReplicates == 1) {
                        myModel.SetSeed(seed);
                    } else {
                        myModel.SetSeed(replicateId, seed);
                    }
                } else {
                    myModel.SetSeed(replicateId);
                }

                // Set the logging behaviour
                myModel.SetExtraSimulationInfo(new String[]{"ReplicateId", "InitSize", "RFrac"},
                        new double[]{replicateId, initialSizeProp, rFrac});
                outFName = outDir + "RepId_" + replicateId + ".csv";
                if (imageFrequency > 0) {
                    myModel.visualiseB = true;
                    String currImgOutDir = imageOutDir + "RepId_" + replicateId;
                    new File(currImgOutDir).mkdirs();
                    myModel.ConfigureImaging(currImgOutDir, imageFrequency);
                }

                // Initialise the simulation
                myModel.InitialiseCellLog(outFName, dt);
                myModel.SetInitialState(initPopSizeArr);
            } else {
                myModel = LoadState(currSavedModelFileName);
            }

            // Run the simulation
            myModel.SetTreatmentSchedule(treatmentScheduleList);
            myModel.Run();
            myModel.Close();
            if (saveModelState) {myModel.SaveModelState(currSavedModelFileName);}
        }
    }

    public static void main(String[] args) {
        int exitCode = new CommandLine(new OnLatticeModel()).execute(args); 
        System.exit(exitCode);
    }
}