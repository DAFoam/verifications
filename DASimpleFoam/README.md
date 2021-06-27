# DASimpleFoam

To compute the adjoint derivatives, run:

<pre>
python runScript.py
</pre>

To compute the reference derivatives using the forward mode AD for a specific design variable, run: 

<pre>
python runScript.py --task=runForwardAD --dvName="shape" --seedIndex=0
</pre>

The above command will run the primal solver with the forward mode AD, and return the derivative for the 0th "shape" design variable. One can follow a similar syntax for other desgin variables and indices.
