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
    public double sp = 1e-3;                        // passenger mutation birth multiplier
    public double sd = 0.1;                         // driver mutation birth multiplier
    public double birth_rate = 0.5;                 // baseline birth rate
    public double death_rate = 0.5;                 // baseline death rate
    public int r0 = 20;                             // initial size of population (r0 x r0)
    public int delete_thresh = 25;                  // ignore clone sizes smaller than this in Muller plots
    public int ExpectedNumberOfClones = 1000000;    // used to create vectors to store clonal information



    // DO NOT CHANGE THESE
    // assumes model is initialized with all cell's progenyID = 1
    public int KpMAX = 0, KdMAX = 1, progenyNextID = 2, sideLen;
    public int[] progenyToParentIDs,driver_status,driver_number,passenger_number,neighborhood=MooreHood(false);
    public int[][] muller;
    Rand rn=new Rand();

    // constructor from parameters
    TumorEvolution(int sideLenth, int totalTime, int modifier){
        super(sideLenth,sideLenth, Cell.class,false,false);

        System.out.println("Building vectors to store clones (this may take a second...)");
        this.progenyToParentIDs = new int[ExpectedNumberOfClones];
        this.driver_status = new int[ExpectedNumberOfClones];
        this.driver_number = new int[ExpectedNumberOfClones];
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
    void Step(){
        for (Cell c:this) {
            c.Birth();
            c.Death();
        }
        CleanShuffle(rn); // shuffle order of agents which birth/death is done each time step
    }

    /*

        main()
            - change the total time, side length, and modifier here

     */

    public static void main(String[] args) {

        int sideLength = 1000;
        int modifier = 50;                  // when to save data (saves when time % modifier == 0)
        int totalTime = 1000;               // total simulation time

        // some folder / file names
        String foldername = "data-output/";
        String baseFilename =  foldername + "tumor_evolution";
        String gifFilename = baseFilename + ".gif";

        TumorEvolution model = new TumorEvolution(sideLength,totalTime,modifier);

        // VISUALIZE
        UIWindow win = new UIWindow("Tumor Evolution", true);
        UIGrid Vis = new UIGrid(sideLength,sideLength, 1);
        GifMaker myGif = new GifMaker(gifFilename, 100,true);
        win.AddCol(0, Vis);
        win.RunGui();

        for (int i = 0; i <= totalTime; i++) {

            model.Step();

            if (i % modifier == 0) {
                DrawCells(model, Vis);
                myGif.AddFrame(Vis);
                System.out.println("time: " + i + " pop: " + model.Pop());

                for (Cell c : model) {
                    model.muller[c.progenyID][(i / modifier)] = model.muller[c.progenyID][(i / modifier)] + 1;
                }
            }
        }


        System.out.println("Reducing phylogenies (this may take a minute)...");
        ReduceParents(foldername, model.muller, model.progenyToParentIDs, model.driver_status, model.ExpectedNumberOfClones, totalTime, modifier, model.delete_thresh);
        System.out.println("Simulation finished...");

        myGif.Close();
        win.Close();

        return;
    }


    /*

        DrawCells()
            - Real Time Visualization of Cells, color-coded by driver number

     */

    public static void DrawCells(TumorEvolution model, UIGrid visCells) {
        for (int i = 0; i < visCells.length; i++) {
            Cell c=model.GetAgent(i);
            visCells.SetPix(i, (c==null) ? RGB(1,1,1) : CategorialColor((c.kd - 1) % 19 ));
        }
    }


    /*

        ReduceParents()
            - removes clones that are dead and gone
            - reconnects parents together so that there is a continuous lineage when parent dies and child lives on
            - writes everything out

     */

    public static void ReduceParents(String foldername, int[][] mullerGenetic, int[] progenyToParentIDs, int[] driver_status, int expectedProgeny, int totalTime, int modifier, int delete_thresh) {

        FileIO parentsReduced = new FileIO((foldername + "parents.csv"), "w");
        FileIO geneticMullerReduced = new FileIO((foldername + "clones.csv"), "w");
        FileIO driverStatusReduced = new FileIO((foldername + "driverStatus.csv"), "w");

        int rowmax = 0;
        ArrayList<Integer> reducedParents = new ArrayList<>();
        ArrayList<Integer> reducedProgeny = new ArrayList<>();
        ArrayList<Integer> reducedDriverStatus = new ArrayList<>();
        reducedParents.add(new Integer(-1)); // add 0th parent
        reducedProgeny.add(new Integer(0)); // add 0th progeny

        for (int jj = 1; jj < expectedProgeny; jj++) {
            rowmax = 0;
            for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                rowmax = (mullerGenetic[jj][ii] > rowmax) ? mullerGenetic[jj][ii]: rowmax;
            }

            if (rowmax > delete_thresh) {

                int mySupposedParent = progenyToParentIDs[jj];

                // check if my parent is in there:
                while (!reducedProgeny.contains(new Integer(mySupposedParent))) {

                    // parent is parent of my parent, lol -- is that parent in there?
                    mySupposedParent = progenyToParentIDs[progenyToParentIDs[mySupposedParent]];

                    // breaks if it is in there, eventually up the chain
                }

                // if my parent is in the reduced thing, add me.
                reducedParents.add(new Integer(reducedProgeny.indexOf(mySupposedParent)));
                reducedProgeny.add(new Integer(jj));
                reducedDriverStatus.add(new Integer(driver_status[jj]));

            }
        }


        // first, output time:
        StringBuilder sb = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        StringBuilder sb3 = new StringBuilder(); //ryan
        for (int time = 0; time < (totalTime / modifier) - 1; time++) { sb.append(time*modifier + ","); }
        sb.append(modifier*((totalTime / modifier) - 1) + "\n");
        geneticMullerReduced.Write(sb.toString()); // file name


        // out gene muller big to file
        int ii_reduced = 0;
        for (int jj = 1; jj < expectedProgeny; jj++) {
            sb = new StringBuilder();

            // check if it's in reduced
            if (reducedProgeny.contains(new Integer(jj))) {
                for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                    sb.append(mullerGenetic[jj][ii] + ",");
                }
                sb.append(mullerGenetic[jj][(totalTime / modifier) - 1] + "\n");

                // output string, including newline \n
                geneticMullerReduced.Write(sb.toString());


                // first iteration no comma
                sb2.append(( (ii_reduced == 0) ? Integer.toString(reducedParents.get(ii_reduced) + 1) : "," + Integer.toString(reducedParents.get(ii_reduced) + 1) )); // add supposed parent, not actual

                // first one is driver, otherwise check driver status
                sb3.append(( (ii_reduced == 0) ? Integer.toString(1) : "," + Integer.toString(reducedDriverStatus.get(ii_reduced)) ));

                ii_reduced++;
            }


        }
        parentsReduced.Write(sb2.toString());
        driverStatusReduced.Write(sb3.toString());
        driverStatusReduced.Close();
        geneticMullerReduced.Close();
        parentsReduced.Close();

    }

}

