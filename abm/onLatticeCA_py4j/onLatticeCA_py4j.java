package onLatticeCA_py4j;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.UIGrid;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;
import py4j.GatewayServer;

import java.util.Arrays;

import static HAL.Util.LoadState;
import static HAL.Util.SaveState;

public class onLatticeCA_py4j extends AgentGrid2D<Cell> implements SerializableModel {
    /* ========================================================================
     * --- Parameters ----
     * ======================================================================== */
    int[] initPopSizeArr; // Initial population sizes
    double divisionRate_S; // Proliferation rate of sensitive cells in d^-1. Proliferation will be attempted at this rate.
    double divisionRate_R; // Proliferation rate of resistant cells in d^-1
    double movementRate_S; // Movement rate of sensitive cells in d^-1. Movement will be attempted at this rate.
    double movementRate_R; // Movement rate of resistant cells in d^-1.
    double deathRate_S; // Natural death rate of sensitive cells in d^-1. Mean life span of a cell will be 1/k
    double deathRate_R;// Natural death rate of resistant cells in d^-1.
    double drugKillProportion; // Drug induced death rate in d^-1.
    double currDrugConcentration; // Current drug concentration in [0,1].
    double[][] treatmentScheduleList; // Treatment schedule in format {{tStart, tEnd, drugConcentration}}
    double tEnd; // End time in d
    double dt; // Time step in d
    int nTSteps; // Number of time steps to simulate
    int[] hood = Util.VonNeumannHood(false); // Define division in von Neumann neighbourhood. Argument indicates that we don't want the center position returned (can only divide into neighbouring squares).
    String initialSeedingType;
    double initRadius = 10.; // Radius circle if cells seeded as circle
    int nCycles = 0; // Number of completed AT cycles; doesn't do anything atm, left in there in case it might be useful in the future.
    int nFailedDivs_R = 0;
    int nAttemptedDivs_R = 0;
    int nDeath_R = 0;

    // Helper variables
    int tIdx = 0;
    int scaleFactor = 2; // Scale factor used to display grid
    int pause = 0; // Pause between time steps to smoothen simulation
    int[] cellCountsArr = new int[2]; // Array to hold number of each cell type
    double[] paramArr;
    Rand rn = new Rand(1); // Random number generator
    Rand rn_ICs = new Rand(1); // Random number generator for initial conditions
    UIGrid vis;
    final static int BLACK = Util.RGB(0, 0, 0);

    // Output
    int verboseLevel = 2;
    Boolean terminateAtProgression = true; // Terminates AT profiling at progression. If not it will be terminated at tEnd
    double printFrequency = 1; // Frequency at which output is printed to the screen
    String cellCountLogFileName;
    FileIO cellCountLogFile = null;
    double logCellCountFrequency = -1; // Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging
    Boolean profilingMode = false; // Whether to log directly to file or save in memory until the end of the simulation
    double[][] outputArr; // Array to hold the simulation results if they're not directly logged to file
    Boolean visualiseB = true;
    int imageFrequency = -1; // Frequency at which an image of the tumour is saved. Negative number turns it off
    private String imageOutDir = "."; // Directory which to save images to
    double[] extraSimulationInfo;
    String[] extraSimulationInfoNames;

    // ------------------------------------------------------------------------------------------------------------
    public static void main(String[] args) {
        // Allow the port to be changed via the command line
        int port = 25333;
        int xDim = 100;
        int yDim = 100;
        for (int i = 0; i < args.length; i+=2) {
            if (args[i].equalsIgnoreCase("-port")) {
                port = Integer.parseInt(args[i + 1]);
            }
            if (args[i].equalsIgnoreCase("-xDim")) {
                xDim = Integer.parseInt(args[i + 1]);
            }
            if (args[i].equalsIgnoreCase("-yDim")) {
                yDim = Integer.parseInt(args[i + 1]);
            }
        }

        // Change to headless mode to allow for faster run time
        System.setProperty("java.awt.headless", "true");

        // Start py4j gateway server
        new GatewayServer(new onLatticeCA_py4j(xDim, yDim), port).start();
        System.out.println("running");
    }

    // ------------------------------------------------------------------------------------------------------------
    // Grid Constructors
    public onLatticeCA_py4j() {
        super(100, 100, Cell.class);
    }

    public onLatticeCA_py4j(int x, int y) {
        super(x, y, Cell.class);
    }

    public onLatticeCA_py4j(int x, int y, double[] paramArr, double dt) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
    }

    public onLatticeCA_py4j(int x, int y, double[] paramArr, double dt, String cellCountLogFileName) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
        this.cellCountLogFileName = cellCountLogFileName;
    }

    public onLatticeCA_py4j(int x, int y, double[] paramArr, double dt, UIGrid vis) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = dt;
        this.vis = vis;
    }

    public onLatticeCA_py4j(int x, int y, double[] paramArr, double dt, UIGrid vis, String cellCountLogFileName) {
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

    public void SetTimeStep(double dt) {this.dt = dt;}

    public void SetInitialState(int[] initialStateArr) {this.initPopSizeArr = initialStateArr;}

    public void SetInitialState(double[] initialStateArr) {
        this.initPopSizeArr = new int[] {(int) initialStateArr[0], (int) initialStateArr[1]};
    }

    public void SetInitialState(int[] initialStateArr, String initialConditionType, int x) {
        this.initPopSizeArr = initialStateArr;
        this.initialSeedingType = initialConditionType;
        if (initialConditionType.equals("circle")) {this.initRadius = x;}
    }

    public void SetInitialState(double[] initialStateArr, String initialConditionType, int x) {
        this.SetInitialState(initialStateArr);
        this.initialSeedingType = initialConditionType;
        if (initialConditionType.equals("circle")) {this.initRadius = x;}
    }

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
                deathRate_S, deathRate_R, drugKillProportion, dt, nCycles, nAttemptedDivs_R, nFailedDivs_R, nDeath_R};
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

    public void InitSimulation_Circle(double radius, double rFrac) {
        //places tumor cells in a circle
        int[] circleHood_all = Util.CircleHood(true, radius);
        int[] circleHood_resistant = Util.CircleHood(true, Math.sqrt(rFrac)*radius); // Inner circle for the resistant cells
        int nCells = MapHood(circleHood_all, xDim / 2, yDim / 2);
        int nCell_R = MapHood(circleHood_resistant, xDim / 2, yDim / 2);
        Arrays.fill(cellCountsArr, 0); // clear the cell counter

        // Fill the inner circle with resistant cells
        for (int i = 0; i < nCell_R; i++) {
            Cell c = NewAgentSQ(circleHood_resistant[i]);
            c.resistance = 1;
            c.divisionRate = this.divisionRate_R;
            c.movementRate = this.movementRate_R;
            c.deathRate = this.deathRate_R;
            c.Draw();
            cellCountsArr[c.resistance] += 1; // Update population counter
        }

        // Fill the outer circle with sensitive cells
        // This is a bit dirty but works. Basically try to place a sensitive cell at all locations
        // within the big circle. If a spot is already taken because there's a resistant cell there
        // then just skip this iteration.
        for (int i = 0; i < nCells; i++) {
            try {
                Cell c = NewAgentSQ(circleHood_all[i]);
                c.resistance = 0;
                c.divisionRate = this.divisionRate_S;
                c.movementRate = this.movementRate_S;
                c.deathRate = this.deathRate_S;
                c.Draw();
                cellCountsArr[c.resistance] += 1;
            } catch (RuntimeException e) {}
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
        nFailedDivs_R = 0;
        nAttemptedDivs_R = 0;
        nDeath_R = 0;
        for (Cell c : this) { //iterate over all cells in the grid
            // Compute the total propensity function (the probability of an event happening)
            totPropensity = (c.divisionRate + c.movementRate +
                    c.deathRate) * dt;

            // Check if an event occurred
            if (rn.Double() < totPropensity) {
                r = rn.Double();
                // 1. Division
                if (r < ((c.divisionRate * dt) / totPropensity)) {
                    nAttemptedDivs_R += c.resistance;
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
                    } else {
                        nFailedDivs_R += c.resistance;}
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
                    nDeath_R += c.resistance;
                }
            }
            cellCountsArr[c.resistance] += 1; // Update population counter
        }
        ShuffleAgents(rn); //shuffle order of for loop iteration over cells
    }

    // ------------------------------------------------------------------------------------------------------------
    public void Run() {
        // Initialise visualisation window
        UIGrid currVis = new UIGrid(xDim,yDim,scaleFactor,visualiseB);
        this.vis = currVis;
        Boolean completedSimulationB = false;
        Boolean logged = false;
        currDrugConcentration = treatmentScheduleList[0][2];
        // Set up the grid and initialise log if this is the beginning of the simulation
        if (tIdx==0) {
            if (initialSeedingType.equalsIgnoreCase("random")) {
                InitSimulation_Random(initPopSizeArr[0], initPopSizeArr[1]);
            } else if (initialSeedingType.equalsIgnoreCase("circle")) {
                InitSimulation_Circle(initRadius,((double) initPopSizeArr[1])/(initPopSizeArr[0]+ initPopSizeArr[1]));
            }
            PrintStatus(0);
            if (cellCountLogFile==null && cellCountLogFileName!=null) {InitialiseCellLog(this.cellCountLogFileName);}
            SaveCurrentCellCount(0);
            SaveTumourImage(tIdx);
            tIdx = 1;
        } else {
            // Continue from a restart
            if (cellCountLogFile==null && cellCountLogFileName!=null) {
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
//        this.Close(logged); xxx
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
            this.profilingMode = profilingMode;
            double[] tmpArr = GetModelState();
            int extraFields = extraSimulationInfoNames==null? 0: extraSimulationInfoNames.length;
            this.outputArr = new double[5][tmpArr.length+extraFields];
            // Initialise the logging array
            for (int i=0; i<outputArr.length; i++) {for (int j=0; j<outputArr[0].length; j++) {outputArr[i][j] = 0;}}
        }
    }

    private void WriteLogFileHeader() {
        cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells,DrugConcentration,rS,rR,mS,mR,dS,dR,dD,dt,NCycles,NAttemptedDivs,NFailedDivs,NDeaths");
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
            if (!profilingMode) {
                cellCountLogFile.WriteDelimit(GetModelState(),",");
                if (extraSimulationInfoNames!=null) {
                    cellCountLogFile.Write(",");
                    cellCountLogFile.WriteDelimit(extraSimulationInfo, ",");
                }
                cellCountLogFile.Write("\n");
                successfulLog = true;
            } else {
                // Make space for new entry in data array
                for (int j=0; j<outputArr[0].length; j++) {
                    for (int i = outputArr.length-1; i > 0; i--) {
                        outputArr[i][j] = outputArr[i-1][j];
                    }
                }

                // Log the new entry
                double[] currModelState = GetModelState();
                for (int j=0; j<currModelState.length; j++) {
                    outputArr[0][j] = currModelState[j];
                }
                if (extraSimulationInfoNames!=null) {
                    for (int j=0; j<extraSimulationInfo.length; j++) {
                        outputArr[0][currModelState.length+j] = extraSimulationInfo[j];
                    }
                }
                successfulLog = true;
            }
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
    public void Reset(Boolean logged) {
        this.Close(logged);
        this.tIdx = 0;
        this.ResetHard();
    }

    public void Close() {
        if (cellCountLogFile!=null) {
            cellCountLogFile.Close();
            cellCountLogFile = null;
        }
    }

    public void Close(Boolean logged) {
        if (!logged) {
            tIdx--;
            SaveCurrentCellCount(0);}
        if (profilingMode) {
            for (int i = 0; i < outputArr.length; i++) {
                cellCountLogFile.WriteDelimit(outputArr[i],",");
                cellCountLogFile.Write("\n");
            }
        }
        this.Close();
    }
}