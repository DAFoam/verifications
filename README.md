# verifications

This repo contains configuration files to verify the accuracy of DAFoam adjoint derivative computation.

Before running, we need to generate the mesh by running:

<pre>
./preProcessing.sh
</pre>

To compute the adjoint derivatives, run:

<pre>
mpirun -np 4 python runScript.py
</pre>

To compute the reference derivatives using the forward mode AD for a specific design variable, run: 

<pre>
mpirun -np 4 python runScript.py --task=runForwardAD --dvName="shape" --seedIndex=0
</pre>

The above command will run the primal solver with the forward mode AD, and print out the derivative for the 0th "shape" design variable to the screen at the end of the computation. One can follow a similar syntax for other design variables and indices.

**NOTE**: If running in serial, (mpirun -np 1), the adjoint derivative is expected to match the forward mode derivative (reference) by about eight significant digits. This is because the default tolerance for flow and adjoint computation is 1e-10. To get more digit match, reduce the above tolerances.

**NOTE**: If running in parallel, the adjoint matches the reference by only six digits, this is because the meshWave function has not been differentiated in parallel. This will be fixed soon.

