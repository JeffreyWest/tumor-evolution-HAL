package Model;

import Framework.GridsAndAgents.AgentGrid2D;
import Framework.Gui.*;
import Framework.Rand;
import static Framework.Util.*;

public class TumorEvolution extends AgentGrid2D<Cell> {

    // parameters for the simulation
    public double mu_p = 1e-4;                      // passenger mutation rate
    public double mu_d = 1e-3;                      // driver mutation rate
    public double sd = 0.1;                         // driver mutation birth multiplier
    public double birth_rate = 0.1;                 // baseline birth rate
    public double death_rate = 0.09;                // baseline death rate
    public int r0 = 20;                             // initial size of population (r0 x r0)
    public int delete_threshold = 25;               // ignore clone sizes smaller than this in Muller plots

    // utility parameters/functions
    public int[] neighborhood=MooreHood(false);
    Rand rn=new Rand();
    Clone clone0;

    // constructor
    TumorEvolution(int sideLenth){
        super(sideLenth,sideLenth, Cell.class,false,false);

        // "seed" clone for common ancestry (non-existent in simulation)
        clone0 = new Clone(null,1,0);

        // construct original tumor
        for (int x = 0; x < r0; x++) {
            for (int y = 0; y < r0; y++) {

                // initialize all new cells w/ 1 driver & 0 passengers
                Clone newClone = new Clone(clone0,1,0);
                NewAgentSQ(x + xDim/2 - r0 / 2,y + yDim/2 - r0 / 2).Init(newClone);
            }
        }
    }

    // Step function for model ("steps" all cells through birth/death/mutation)
    void Step(){

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

        int sideLength = 50;    // domain size
        int modifier = 100;     // when to save data (saves when time % modifier == 0)
        int totalTime = 5000;   // total simulation time

        // some folder / file names
        String foldername = "data-output/";
        TumorEvolution tumor = new TumorEvolution(sideLength);
        GridWindow window = new GridWindow("Tumor Evolution", sideLength,sideLength,5,true);
        GifMaker myGif = new GifMaker(foldername + "tumor_evolution.gif", 250,true);

        for (int time = 0; time <= totalTime; time++) {

            tumor.Step(); // birth & death for all cells

            if (time % modifier == 0) {
                DrawCellsAndSaveGif(tumor, window, myGif);

                System.out.println(time);
                // record clonal pops
                tumor.clone0.RecordClones(time);  // iterate through all clones & record population sizes

            }
        }

        // save clonal information in EvoFreq format:
        String[] AttributesList = new String[]{"Drivers","Passengers,Color"};
        tumor.clone0.PopRecordToCSV(foldername+"phylogeny_tracker.csv",AttributesList, (Clone c) -> {return GetAttributes(c); },tumor.delete_threshold);

        myGif.Close();
        window.Close();

        return;
    }


    /*

        DrawCells()
            - Real Time Visualization of Cells, color-coded by driver number

     */

    public static void DrawCellsAndSaveGif(TumorEvolution model, GridWindow vis, GifMaker myGif) {
        for (int i = 0; i < vis.length; i++) {
            Cell c=model.GetAgent(i);
            vis.SetPix(i, (c==null) ? RGB(1,1,1) : c.clone.color);
        }

        myGif.AddFrame(vis);
    }

    /*

        GetAttributes()
            - get clone-specific attributes to output to EvoFreq csv

     */

    public static String[] GetAttributes(Clone c) {
        return new String[]{Integer.toString(c.kd),Integer.toString(c.kp),c.Hex()};
    }

}

