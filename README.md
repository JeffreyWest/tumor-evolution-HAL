## Tissue structure accelerates evolution: premalignant sweeps precede neutral invasion

This is the open-source code repository for the passenger driver model of tumor evolution. We encourage you to begin by reading the manuscript here: (link coming soon). There is a two-dimensional on lattice agent-based model called PD2D and a three-dimensional on lattice agent-based model of tumor evolution for ductal carcinoma in situ.

## Before you start
TumorEvolution.java uses the open-source agent-based modeling platform called [HAL](http://halloworld.org/). In order to run this code base, you'll need to download the latest version of [Java](http://www.oracle.com/technetwork/java/javase/downloads/jdk9-downloads-3848520.html) and an editor (we suggest using [IntelliJ Idea](https://www.jetbrains.com/idea/download/)).

### Setting up the project in IntelliJ Idea

1. Open Intellij Idea and click "create project from existing sources" ("file/ new/ project from existing sources" from the main GUI) and direct it to the unzipped AgentFramework Source code directory.
2. Continue through the rest of the setup, click next until it asks for the Java SDK:
- "/Library/ Java/ JavaVirtualMachines/" on Mac.
- "C:\ Program Files\ Java\" on Windows.
3. Once the setup is complete we will need to do one more step and add some libraries that allow for 2D and 3D OpenGL visualization:
4. open the Intellij IDEA main gui
5. go to "file/ project structure"
6. click the "libraries" tab
7. use the minus button to remove any pre-existing library setup
8. click the plus button, and direct the file browser to the "Framework/ lib" folder.
9. click apply or ok
