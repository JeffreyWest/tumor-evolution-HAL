package Model;

import Framework.GridsAndAgents.AgentGrid2D;
import Framework.Gui.*;
import Framework.Tools.FileIO;
import Framework.Rand;
import static Framework.Util.*;
import java.util.ArrayList;

public class TumorEvolution extends AgentGrid2D<Cell> {

    // parameters for the simulation
    public double mu_p = 1e-3;                      // passenger mutation rate
    public double mu_d = 1e-4;                      // driver mutation rate
    public double sd = 0.1;                         // driver mutation birth multiplier
    public double birth_rate = 0.5;                 // baseline birth rate
    public double death_rate = 0.5;                 // baseline death rate
    public int r0 = 20;                             // initial size of population (r0 x r0)
    public int delete_threshold = 25;               // ignore clone sizes smaller than this in Muller plots
    public int ExpectedNumberOfClones = 1000000;    // used to create vectors to store clonal information

    // DO NOT CHANGE THESE (assumes model is initialized with all cell's progenyID = 1)
    public int KpMAX = 0, KdMAX = 1, progenyNextID = 2, sideLen;
    public int[] progenyToParentIDs,driver_status,passenger_number,neighborhood=MooreHood(false);
    public int[][] muller;
    Rand rn=new Rand();

    // constructor from parameters
    TumorEvolution(int sideLenth, int totalTime, int modifier){
        super(sideLenth,sideLenth, Cell.class,false,false);

        System.out.println("Building vectors to store clones (this may take a second...)");
        this.progenyToParentIDs = new int[ExpectedNumberOfClones];
        this.driver_status = new int[ExpectedNumberOfClones];
        this.passenger_number = new int[ExpectedNumberOfClones];
        this.sideLen = sideLenth;

        // construct muller vector:
        muller = new int[ExpectedNumberOfClones][totalTime / modifier];
        for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {
            for (int iii = 0; iii < ExpectedNumberOfClones; iii++) {
                muller[iii][jjj] = 0;
            }
        }

        // construct original tumor
        for (int x = 0; x < r0; x++) {
            for (int y = 0; y < r0; y++) {
                NewAgentSQ(x + sideLen/2 - r0 / 2,y + sideLen/2 - r0 / 2).Init(0,1,1,0);
            }
        }
    }

    // Step function for model ("steps" all cells through birth/death/mutation)
    void Step(int time, int modifier){

        // store muller information, if the time is right
        if (time % modifier == 0) {
            for (Cell c : this) {
                this.muller[c.progenyID][(time / modifier)] = this.muller[c.progenyID][(time / modifier)] + 1;
            }
        }

        // step through birth/death process for all cells
        for (Cell c:this) {
            c.Birth();
            c.Death();
        }

        // shuffle order of agents which birth/death is done each time step
        CleanShuffle(rn);
    }

    /*

        main()
            - change the total time, side length, and modifier here

     */

    public static void main(String[] args) {

        int sideLength = 500;               // domain size
        int modifier = 50;                  // when to save data (saves when time % modifier == 0)
        int totalTime = 500;                // total simulation time

        // some folder / file names
        String foldername = "data-output/";
        String baseFilename =  foldername + "tumor_evolution";
        String gifFilename = baseFilename + ".gif";

        TumorEvolution model = new TumorEvolution(sideLength,totalTime,modifier);

        // VISUALIZE
        UIWindow win = new UIWindow("Tumor Evolution", true);
        UIGrid Vis = new UIGrid(sideLength,sideLength, 1);
        win.AddCol(0, Vis);
        win.RunGui();

        GifMaker myGif = new GifMaker(gifFilename, 100,true);

        for (int time = 0; time < totalTime; time++) {

            model.Step(time,modifier);

            if (time % modifier == 0) {
                DrawCellsAndSaveGif(model, Vis, myGif);
                System.out.println("time: " + time + " pop: " + model.Pop());
            }
        }

        System.out.println("Reducing phylogenies (this may take a minute)...");
        model.ReduceAndOutputClones(foldername, totalTime, modifier);
        System.out.println("Simulation finished...");

        myGif.Close();
        win.Close();

        return;
    }


    /*

        DrawCells()
            - Real Time Visualization of Cells, color-coded by driver number

     */

    public static void DrawCellsAndSaveGif(TumorEvolution model, UIGrid vis, GifMaker myGif) {
        for (int i = 0; i < vis.length; i++) {
            Cell c=model.GetAgent(i);
            vis.SetPix(i, (c==null) ? RGB(1,1,1) : CategorialColor((c.kd - 1) % 19 ));
        }

        myGif.AddFrame(vis);
    }


    /*

        ReduceParents()
            - removes clones that never grew above "delete_threshold"
            - reconnects parents together so that there is a continuous lineage when parent dies and child lives on
            - writes everything out

     */

    public void ReduceAndOutputClones(String foldername, int totalTime, int modifier) {

        FileIO mullerReduced = new FileIO((foldername + "evofreq_dataframe.csv"), "w");

        ArrayList<Integer> reducedParents = new ArrayList<>();
        ArrayList<Integer> reducedProgeny = new ArrayList<>();
        ArrayList<Integer> reducedDriverStatus = new ArrayList<>();
        reducedParents.add(new Integer(-1)); // add 0th parent
        reducedProgeny.add(new Integer(0)); // add 0th progeny

        for (int jj = 1; jj < this.ExpectedNumberOfClones; jj++) {
            int rowmax = 0;
            for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                rowmax = (muller[jj][ii] > rowmax) ? muller[jj][ii]: rowmax;
            }

            if (rowmax > delete_threshold) {

                int mySupposedParent = progenyToParentIDs[jj];

                // check if my parent is in there:
                // (it may not be, if parent is below delete_threshold)
                while (!reducedProgeny.contains(new Integer(mySupposedParent))) {

                    // parent is parent of my parent, etc
                    mySupposedParent = progenyToParentIDs[progenyToParentIDs[mySupposedParent]];

                    // loop breaks when a parent is found eventually up the chain of parents
                }

                // add the (parent, progeny) pair to respective array lists
                reducedParents.add(new Integer(reducedProgeny.indexOf(mySupposedParent)));
                reducedProgeny.add(new Integer(jj));

                // also add the driver status to its array list
                reducedDriverStatus.add(new Integer(driver_status[jj]));

            }
        }

        StringBuilder sb = new StringBuilder();

        // add a row of headers (including each time point) to the first row of .csv file
        sb.append("CloneID,Parent,Drivers,");
        for (int time = 0; time < (totalTime / modifier) - 1; time++) { sb.append(time*modifier + ","); }
        sb.append(modifier*((totalTime / modifier) - 1) + "\n");
        mullerReduced.Write(sb.toString()); // file name

        int ii_reduced = 0;
        for (int jj = 1; jj < this.ExpectedNumberOfClones; jj++) {
            sb = new StringBuilder();

            // check if it's in reduced array list (i.e. that it exceeded delete_threshold)
            if (reducedProgeny.contains(new Integer(jj))) {

                // append cloneID, parent and # of drivers:
                sb.append(Integer.toString(ii_reduced)+ "," + Integer.toString(reducedParents.get(ii_reduced) + 1)+ ",");
                sb.append(( (ii_reduced == 0) ? Integer.toString(1) : Integer.toString(reducedDriverStatus.get(ii_reduced)) ));
                sb.append(",");

                // append this cloneID's size over time
                for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                    sb.append(muller[jj][ii] + ",");
                }
                sb.append(muller[jj][(totalTime / modifier) - 1] + "\n");

                // write out
                mullerReduced.Write(sb.toString());
                ii_reduced++;
            }
        }

        // close file
        mullerReduced.Close();
    }
}

