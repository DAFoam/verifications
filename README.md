# verifications

This repo contains configuration files to verify the accuracy of DAFoam adjoint derivative computation (DAFoam achieves **machine precision accurate adjoint**). Here we compared the derivatives between the forward mode AD (reference)and the JacobianFree adjoint method.

Before running, we need to generate the mesh by running:

<pre>
./preProcessing.sh
</pre>

To compute the adjoint derivatives, run:

<pre>
mpirun -np 4 python runScript.py --mode=reverse --task=runAdjoint
</pre>

To compute the reference derivatives using the forward mode AD for a specific design variable, run: 

<pre>
mpirun -np 4 python runScript.py --mode=forward --task=runForwardAD --dvName="shape" --seedIndex=0
</pre>

The above command will run the primal solver with the forward mode AD, and print out the derivative for the 0th "shape" design variable to the screen during the computation. One can follow a similar syntax for other design variables and indices.

