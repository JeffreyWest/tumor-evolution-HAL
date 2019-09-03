package Framework.Tools.PhylogenyTracker;

import Framework.Interfaces.GetGenomeAttrs;
import Framework.Tools.FileIO;
import Framework.Util;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Created by bravorr on 8/4/17.
 */
class PhylogenyTrackerInternal<T extends Genome> implements Iterable<T> {
    int nGenomesEver;
    int nLivingGenomes;
    int nTreeGenomes;
    long totalPop;
    final public boolean removeEmptyLeaves;
    T progenitor;
    T listFirst;
    ArrayList<long[]>popRecords;
    ArrayList<String> recordLabels;

    public PhylogenyTrackerInternal(T progenitor, boolean removeEmptyLeaves) {
        this.progenitor = progenitor;
        this.removeEmptyLeaves = removeEmptyLeaves;
        this.listFirst=progenitor;
        ResetRecord();
    }

    public void ResetRecord(){
        popRecords=new ArrayList<>();
        recordLabels=new ArrayList<>();
    }

    @Override
    public Iterator<T> iterator() {
        return new myIter(progenitor);
    }

    private class myIter implements Iterator<T> {
        T curr;

        myIter(T last) {
            this.curr = last;
        }

        @Override
        public boolean hasNext() {
            return curr != null;
        }

        @Override
        public T next() {
            T ret=curr;
            curr = (T) curr.prev;
            return ret;
        }
    }
    public void RecordClones(int timepoint){
        RecordClones(Integer.toString(timepoint));
    }

    public void RecordClones(String timepointLabel){
        if(removeEmptyLeaves){
            throw new IllegalStateException("can't record pops with leaf removal on!");
        }
        if(timepointLabel=="CloneID"||timepointLabel=="ParentID"){
            throw new IllegalArgumentException("label cannot be CloneID or ParentID!");
        }
        long[]pops=new long[nGenomesEver];
        int ct=0;
        for (Object cloneobj : progenitor) {
            Genome clone=(Genome)cloneobj;
            pops[clone.id]=clone.GetPop();
        }
        popRecords.add(pops);
        recordLabels.add(timepointLabel);
    }

    public void PopRecordToCSV(String path, String[]AttrHeaders, GetGenomeAttrs<T> GetAttrs,int includePopCutoff){
        FileIO out=new FileIO(path,"w");
        for (String header : AttrHeaders) {
            if(header=="CloneID"||header=="ParentID"){
                throw new IllegalArgumentException("attr header cannot be CloneID or ParentID!");
            }
        }
        out.Write(Util.ArrToString(AttrHeaders,",")+",CloneID,ParentID,"+Util.ArrToString(recordLabels,",")+"\n");
        long[]maxpops=new long[nGenomesEver];
        int[]parentIDs=new int[nGenomesEver];
        String[]attrsOut=new String[nGenomesEver];
        int id=0;
        for (Object cloneobj : progenitor) {
            T clone=(T)cloneobj;
            int parentid=0;
            if(id!=0){
                parentid=clone.GetParent().id;
            }
            parentIDs[id]=parentid;
            attrsOut[id]=Util.ArrToString(GetAttrs.GetAttrs(clone),",");
            id++;
        }
        //get maxpops
        for (long[] pops : popRecords) {
            for (int i = 0; i < pops.length; i++) {
                maxpops[i]=Math.max(maxpops[i],pops[i]);
            }
        }
        //propagate maxpops up to see what nodes will stay
        for(int i=maxpops.length-1;i>=0;i--){
            maxpops[parentIDs[i]]=Math.max(maxpops[parentIDs[i]],maxpops[i]);
        }
        for (id = 0; id <nGenomesEver ; id++) {
            if(maxpops[id]>includePopCutoff) {
                out.Write(attrsOut[id] + "," + id + "," + parentIDs[id] + ",");
                //get pops at each timepoint
                for (int j = 0; j < popRecords.size(); j++) {
                    long pop = 0;
                    if (popRecords.get(j).length > id) {
                        pop = popRecords.get(j)[id];
                    }
                    if (j < popRecords.size() - 1) {
                        out.Write(pop + ",");
                    } else {
                        out.Write(pop + "");
                    }
                }
                out.Write("\n");
            }
        }
        out.Close();
    }

    public int GetNumGenomes() {
        return nGenomesEver;
    }

    public int GetNumLivingGenomes() {
        return nLivingGenomes;
    }

    public int GetNumTreeGenomes() {
        return nTreeGenomes;
    }

    public long GetTotalPop() {
        return totalPop;
    }

    public T GetProgentior(){
        return progenitor;
    }
}
