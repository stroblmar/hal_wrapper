package onLatticeCA_jar;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.UIGrid;
import HAL.Gui.GridWindow;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;

import java.io.File;
import java.util.Arrays;

import static HAL.Util.LoadState;
import static HAL.Util.SaveState;

public class OnLatticeGrid extends AgentGrid2D<Cell> implements SerializableModel {
    /* ========================================================================
     * --- Parameters ----
     * ======================================================================== */
    int xDim = 100; // x-dimension of domain (in lattice sites)
    int yDim = 100; // y-dimension of domain (in lattice sites)
    int[] initPopSizeArr; // Initial population sizes
    double initialSize = 1000; // Total initial cell number (absolute)
    double initialSizeProp; // Initial cell density relative to (physical) carrying capacity
    double rFrac = 0.05; // Initial resistance fraction in [0,1]
    double divisionRate_S = 0.027; // Proliferation rate of sensitive cells in d^-1. Proliferation will be attempted at this rate.
    double divisionRate_R = 0.027; // Proliferation rate of resistant cells in d^-1
    double movementRate_S = 0; // Movement rate of sensitive cells in d^-1. Movement will be attempted at this rate.
    double movementRate_R = 0; // Movement rate of resistant cells in d^-1.
    double deathRate_S = 0.0027; // Natural death rate of sensitive cells in d^-1. Mean life span of a cell will be 1/k
    double deathRate_R = 0.0027;// Natural death rate of resistant cells in d^-1.
    double drugKillProportion = 0.75; // Drug induced death rate in d^-1.
    int[] hood = Util.VonNeumannHood(false); // Define division in von Neumann neighbourhood. Argument indicates that we don't want the center position returned (can only divide into neighbouring squares).
    double tEnd = 7; // End time in d
    double dt = 1; // Time step in d
    int nTSteps; // Number of time steps to simulate
    double[][] treatmentScheduleList = {{0, tEnd, 0}}; // Treatment schedule in format {{tStart, tEnd, drugConcentration}}
    int nReplicates = 10; // Number of replicates
    
    // Helper variables
    int tIdx = 0;
    double currDrugConcentration; // Current drug concentration in [0,1].
    int[] cellCountsArr = new int[2]; // Array to hold number of each cell type
    double[] paramArr;
    int seed = 0; // Random number seed
    Rand rn; // Random number generator
    Rand rn_ICs; // Random number generator for initial conditions
    
    // Output - Text
    int verboseLevel = 2;
    double printFrequency = 1; // Frequency at which output is printed to the screen
    String cellCountLogFileName;
    FileIO cellCountLogFile = null;
    double[][] outputArr; // Array to hold the simulation results if they're not directly logged to file
    double logCellCountFrequency = -1; // Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging
    double[] extraSimulationInfo;
    String[] extraSimulationInfoNames;

    // Output - Visualisation
    UIGrid vis;
    Boolean visualiseB = true; // Whether or not to show visualization
    int scaleFactor = 2; // Scale factor used to display grid
    int pause = 0; // Pause betwween time steps to smoothen simulation
    final static int BLACK = Util.RGB(0, 0, 0);
    int imageFrequency = -1; // Frequency at which an image of the tumour is saved. Negative number turns it off
    String imageOutDir = ""; // Directory which to save images to

    // Output - Model file
    Boolean saveModelState = false;
    Boolean fromScratch = true;
    String savedModelFileName = null; // Name of model file to load when continuing a previous run
    // ------------------------------------------------------------------------------------------------------------
    public static void main(String[] args) {
        // Runs an example simulation using the default parameters
        // Setup
        OnLatticeGrid myModel = new OnLatticeGrid();
        int nReplicates = 1; 
        int[] replicateIdList = new int[] {myModel.seed};
        int replicateId;
        int[] initPopSizeArr = {0, 0};
        double[][] treatmentScheduleList = {{0, 50, 0}}; // Treatment schedule in format {{tStart, tEnd, drugConcentration}}
        String outDir = "./tmp/";
        Boolean fromScratch = true; 
        Boolean visualiseB = true;
        int pause = 250;
        String savedModelFileName = null; // Name of model file to load when continuing a previous run
        String currSavedModelFileName;

        // Update model object with new params
        myModel.fromScratch = fromScratch;
        
        // Run the sweep
        // System.setProperty("java.awt.headless", "true"); xxx
        // Organise replicates. If we only do one, it will have the id equals to the random number seed.
        // Otherwise, they'll be numbered 0 to nReplicates.
        if (nReplicates>1) {
            replicateIdList = new int[nReplicates];
            for (int i=0; i<nReplicates; i++) {replicateIdList[i]=i;}
        }
        if (myModel.initialSizeProp>0) { // Make it so that can use proportional definition of initial density
            myModel.initialSize = (int) (myModel.initialSizeProp * myModel.xDim * myModel.yDim);
        }
        initPopSizeArr = new int[]{(int) Math.floor(myModel.initialSize * (1 - myModel.rFrac)), (int) Math.ceil(myModel.initialSize * myModel.rFrac)};
        new File(outDir).mkdirs();
        for (int replicateIdx = 0; replicateIdx < nReplicates; replicateIdx++) {
            replicateId = replicateIdList[replicateIdx];
            currSavedModelFileName = savedModelFileName==null? outDir + "RepId_" + replicateId + ".bin": savedModelFileName;

            // Models can be created either from scratch, or be loaded from a previous run
            if (fromScratch) {
                // Set the random number seed
                if (nReplicates == 1) {
                    myModel.SetSeed(myModel.seed);
                } else {
                    myModel.SetSeed(replicateId, myModel.seed);
                }

                // Set the logging behaviour
                myModel.ConfigureVisualisation(visualiseB, pause);
                if (myModel.imageOutDir.length() > 0) {
                    myModel.visualiseB = true;
                    String currImgOutDir = outDir + "/img/" + "RepId_" + replicateId;
                    new File(currImgOutDir).mkdirs();
                    myModel.ConfigureImaging(currImgOutDir, myModel.imageFrequency);
                }

                // Initialise the simulation
                myModel.InitialiseCellLog(outDir + "RepId_" + replicateId + ".csv");
                myModel.SetInitialState(initPopSizeArr);
            } else {
                myModel = LoadState(currSavedModelFileName);
            }

            // Run the simulation
            myModel.SetTreatmentSchedule(treatmentScheduleList);
            myModel.Run();
            myModel.Close();
            if (myModel.saveModelState) {myModel.SaveModelState(currSavedModelFileName);}
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    // Grid Constructors
    public OnLatticeGrid() {
        super(100, 100, Cell.class);
    }

    public OnLatticeGrid(int x, int y, double[] paramArr, double dt) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
    }

    public OnLatticeGrid(int x, int y, double[] paramArr, double dt, String cellCountLogFileName) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
        this.cellCountLogFileName = cellCountLogFileName;
    }

    public OnLatticeGrid(int x, int y, double[] paramArr, double dt, UIGrid vis) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
        this.vis = vis;
    }

    public OnLatticeGrid(int x, int y, double[] paramArr, double dt, UIGrid vis, String cellCountLogFileName) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
        this.vis = vis;
        this.cellCountLogFileName = cellCountLogFileName;
    }

    // Function used as part of SerializableModel to allow saving the model's state so that I can restart
    // the simulation exactly where I left off at the end of the last treatment interval.
    @Override
    public void SetupConstructors() {
        _PassAgentConstructor(Cell.class);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Functions to parameterise the model
    public void SetParameters(double[] paramArr) {
        this.divisionRate_S = paramArr[0];
        this.divisionRate_R = paramArr[1];
        this.movementRate_S = paramArr[2];
        this.movementRate_R = paramArr[3];
        this.deathRate_S = paramArr[4];
        this.deathRate_R = paramArr[5];
        this.drugKillProportion = paramArr[6];
    }

    public void SetInitialState(int[] initialStateArr) {this.initPopSizeArr = initialStateArr;}

    public void SetDrugConcentration(double drugConcentration) {
        this.currDrugConcentration = drugConcentration;
    }

    public void SetTreatmentSchedule(double[][] treatmentScheduleList) {
        this.treatmentScheduleList = treatmentScheduleList;
    }

    public void SetSeed(int seed) {
        this.rn = new Rand(seed);
        this.rn_ICs = new Rand(seed);
    }

    public void SetSeed(int seed_Simulation, int seed_ICs) {
        this.rn = new Rand(seed_Simulation);
        this.rn_ICs = new Rand(seed_ICs);
    }

    public void SetVerboseness(int verboseLevel, double printFrequency) {
        this.verboseLevel = verboseLevel;
        this.printFrequency = printFrequency;
    }

    public void ConfigureVisualisation(boolean visualiseB, int pause) {
        this.visualiseB = visualiseB;
        this.pause = pause;
    }

    public void ConfigureImaging(String imageOutDir, int imageFrequency) {
        /*
        * Configure location of and frequency at which tumour is imaged.
        */
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
    }

    public double[] GetParameters() {
        return new double[] {divisionRate_S, divisionRate_R, movementRate_S, movementRate_R,
                deathRate_S, deathRate_R, drugKillProportion};
    }

    public double[] GetModelState() {
        return new double[] {tIdx, tIdx*dt, cellCountsArr[0], cellCountsArr[1], Util.ArraySum(cellCountsArr),
                currDrugConcentration, divisionRate_S, divisionRate_R, movementRate_S, movementRate_R,
                deathRate_S, deathRate_R, drugKillProportion, dt};
    }

    // ------------------------------------------------------------------------------------------------------------
    // Seeding functions for the cells
    public void InitSimulation_Random(int s0, int r0) {
        //places tumor cells randomly on the dish
        int[] distributionOfResistanceArr = new int[s0 + r0];
        int[] initialLocArr = new int[xDim * yDim];
        Arrays.fill(cellCountsArr, 0); // clear the cell counter

        // Generate a list of random grid positions
        for (int i = 0; i < xDim * yDim; i++) {
            initialLocArr[i] = i;
        }
        rn_ICs.Shuffle(initialLocArr);

        // Generate a list of random assignment to sensitive or resistance
        for (int i = 0; i < s0; i++) {
            distributionOfResistanceArr[i] = 0;
        }
        for (int i = s0; i < s0 + r0; i++) {
            distributionOfResistanceArr[i] = 1;
        }
        rn_ICs.Shuffle(distributionOfResistanceArr);

        // Create cells according to the above assignments
        for (int i = 0; i < s0 + r0; i++) {
            Cell c = NewAgentSQ(initialLocArr[i]);
            c.resistance = distributionOfResistanceArr[i];
            if (c.resistance == 0) {
                c.divisionRate = this.divisionRate_S;
                c.movementRate = this.movementRate_S;
                c.deathRate = this.deathRate_S;
            } else {
                c.divisionRate = this.divisionRate_R;
                c.movementRate = this.movementRate_R;
                c.deathRate = this.deathRate_R;
            }
            c.Draw();
            cellCountsArr[c.resistance] += 1; // Update population counter
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    public void StepCells() {
        Arrays.fill(cellCountsArr, 0);//clear the cell counts
        double totPropensity;
        int currPos;
        boolean moveSuccessB;
        double deathProb;
        double r;
        ShuffleAgents(rn); //shuffle order of for loop iteration over cells
        
        for (Cell c : this) { //iterate over all cells in the grid
            // Compute the total propensity function (the probability of an event happening)
            totPropensity = (c.divisionRate + c.movementRate +
                    c.deathRate) * dt;

            // Check if an event occured
            if (rn.Double() < totPropensity) {
                r = rn.Double();
                // 1. Division
                if (r < ((c.divisionRate * dt) / totPropensity)) {
                    if (c.HasSpace()) {
                        // If drug is present, check if cell is killed
                        deathProb = drugKillProportion*currDrugConcentration*(1-c.resistance); // Note the kill rate is technically a proportion
                        if (rn.Double() < deathProb) {
                            vis.SetPix(c.Isq(), BLACK); // Death
                            c.Dispose();// Removes cell from spatial grid and iteration
                            cellCountsArr[c.resistance] -= 1; // Update population counter
                        } else {
                            c.Divide();
                            cellCountsArr[c.resistance] += 1;
                        }
                    }
                // 2. Movement
                } else if (r < ((c.divisionRate + c.movementRate) * dt) / totPropensity) {
                    currPos = c.Isq();
                    moveSuccessB = c.Move();
                    if (moveSuccessB) {
                        vis.SetPix(currPos, BLACK); // Remove from old position
                        c.Draw(); // Draw in new place
                    }
                // 3. Natural Death
                } else {
                    vis.SetPix(c.Isq(), BLACK); // Death
                    c.Dispose();// Removes cell from spatial grid and iteration
                    cellCountsArr[c.resistance] -= 1; // Update population counter
                }
            }
            cellCountsArr[c.resistance] += 1; // Update population counter
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    public void Run() {
        // Initialise visualisation window
        // UIGrid currVis = new UIGrid(xDim, yDim, scaleFactor, visualiseB);
        // this.vis = currVis;
        // GridWindow vis=new GridWindow(xDim,yDim,scaleFactor);//used for visualization
        // vis.AddAlphaGrid(currVis);
        UIGrid currVis = new UIGrid(xDim, yDim, scaleFactor, visualiseB);
        this.vis = currVis;
        Boolean completedSimulationB = false;
        Boolean logged = false;
        currDrugConcentration = treatmentScheduleList[0][2];
        // Set up the grid and initialise log if this is the beginning of the simulation
        if (tIdx==0) {
            InitSimulation_Random(initPopSizeArr[0], initPopSizeArr[1]);
            PrintStatus(0);
            if (cellCountLogFile==null && cellCountLogFileName!=null) {InitialiseCellLog(this.cellCountLogFileName);}
            SaveCurrentCellCount(0);
            SaveTumourImage(tIdx);
            tIdx = 1;
        } else {
            // Continue from a restart
            if (cellCountLogFileName!=null) {
                cellCountLogFile = new FileIO(cellCountLogFileName, "a");
            }
            for (Cell c : this) {
                c.Draw();
            }
        }

        // Run the simulation
        double currIntervalEnd;
        int initialCellNumber = Util.ArraySum(cellCountsArr);
        if (treatmentScheduleList==null) treatmentScheduleList = new double[][]{{0,tEnd,currDrugConcentration}};
        for (int intervalIdx=0; intervalIdx<treatmentScheduleList.length;intervalIdx++) {
            currIntervalEnd = treatmentScheduleList[intervalIdx][1];
            nTSteps = (int) Math.ceil(currIntervalEnd/dt);
            currDrugConcentration = treatmentScheduleList[intervalIdx][2];
            completedSimulationB = false;
            while (!completedSimulationB) {
                vis.TickPause(pause);
                StepCells();
                PrintStatus(tIdx);
                logged = SaveCurrentCellCount(tIdx);
                SaveTumourImage(tIdx);
                tIdx++;
                // Check if the stopping condition is met
                completedSimulationB = (tIdx>nTSteps)?true:false;
            }
        }

        // Close the simulation
        this.Close(logged);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Manage and save output
    public void InitialiseCellLog(String cellCountLogFileName) {
        InitialiseCellLog(cellCountLogFileName, 1.);
    }

    public void InitialiseCellLog(String cellCountLogFileName, double frequency) {
        InitialiseCellLog(cellCountLogFileName,frequency,false);
    }

    public void InitialiseCellLog(String cellCountLogFileName, double frequency, Boolean profilingMode) {
        cellCountLogFile = new FileIO(cellCountLogFileName, "w");
        WriteLogFileHeader();
        this.cellCountLogFileName = cellCountLogFileName;
        this.logCellCountFrequency = frequency;
        if (profilingMode) {
            double[] tmpArr = GetModelState();
            int extraFields = extraSimulationInfoNames==null? 0: extraSimulationInfoNames.length;
            this.outputArr = new double[5][tmpArr.length+extraFields];
            // Initialise the logging array
            for (int i=0; i<outputArr.length; i++) {for (int j=0; j<outputArr[0].length; j++) {outputArr[i][j] = 0;}}
        }
    }

    private void WriteLogFileHeader() {
        cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells,DrugConcentration,rS,rR,mS,mR,dS,dR,dD,dt");
        if (extraSimulationInfoNames!=null) {
            cellCountLogFile.Write(",");
            cellCountLogFile.WriteDelimit(extraSimulationInfoNames, ",");
        }
        cellCountLogFile.Write("\n");
    }

    public void SetExtraSimulationInfo(String[] extraInfoNamesArr, double[] extraInfoArr) {
        this.extraSimulationInfoNames = extraInfoNamesArr;
        this.extraSimulationInfo = extraInfoArr;
    }

    public Boolean SaveCurrentCellCount(int currTimeIdx) {
        Boolean successfulLog = false;
        if ((currTimeIdx % (int) (logCellCountFrequency/dt)) == 0 && logCellCountFrequency > 0) {
            cellCountLogFile.WriteDelimit(GetModelState(),",");
            if (extraSimulationInfoNames!=null) {
                cellCountLogFile.Write(",");
                cellCountLogFile.WriteDelimit(extraSimulationInfo, ",");
            }
            cellCountLogFile.Write("\n");
            successfulLog = true;
        }
        return successfulLog;
    }

    public void PrintStatus(int currTimeIdx) {
        if (verboseLevel > 0 && (currTimeIdx % (int) (printFrequency/dt)) == 0) {
            System.out.println("Time: " + currTimeIdx*dt+
                    " - Population Size: " + Util.ArraySum(cellCountsArr)+
                    " - Drug Concentration: " + currDrugConcentration);
        }
    }

    public void SaveTumourImage(int currTimeIdx) {
        if (imageFrequency > 0 && (currTimeIdx % (int) (imageFrequency/dt)) == 0) {
            this.vis.ToPNG(imageOutDir +"/img_t_"+currTimeIdx*dt+".png");
        }
    }

    public void SaveModelState(String stateFileName) {
        // Can't have active pointers when saving the model. So, close everything here.
        if (cellCountLogFile!=null) {
            cellCountLogFile.Close();
            cellCountLogFile = null;
        }
        if (vis!=null) {vis = null;}
        SaveState(this,stateFileName);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Function to clean up loose ends after running a simulation.
    public void Close() {
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }

    public void Close(Boolean logged) {
        if (!logged) {
            tIdx--;
            SaveCurrentCellCount(0);}
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }
}