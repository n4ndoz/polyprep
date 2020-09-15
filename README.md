# polyprep <br>
###########################################<br>
# PolyProtein REpresentation Preparator   #<br>
###########################################<br>
<br>
<br>
## Description<br>
This script is one of the results of my undergrad Work as a research intern
at Inmetro (RJ, Brazil). At the time my work was to study substract recognition
on the HIV-1 protease, in light of machine learning methods. This code was
a supporting partof both my work and my Tutor's (MFR Dias) PhD Thesis.
This code consists of a routine of fragmenting a polyprotein aminoacid 
sequence intelligently, and preparing it for docking experiments (using
Autoodck 4 or Vina). It performs fragmentation first by searching for
given interfaces (on a sepparate fasta file) inside the sequence, and then
gathering fragments of these interfaces of also given frag_lengths. It also
saves the fragments that do not contain interfaces, labeling them as NEG(ative).
If desired, the user can then model these fragments using any template wanted*,
(more on this later) by calling modeller succesively.<br>
<br>
The ideia behind looking for the interfaces on the input sequence(s) and
not only by combining the letters, is to apply this method on sequence space. 
For instance imagine we have frag_1 (seq = APRL*CKVI), where * is the cleavage
interface. If we have a positive interface with seq = RL*CK, passing a sliding
window and gathering all possible 4aa long residues yields:<br>
PRL*C<br>
RL*CK<br>
L*CKI<br>
<br>
This was done because at the time we were dealing with a heavely unballanced 
data set/sequence space. For each len_frag value, we had several folds more 
negative examples than positives, something that was cripling attempts on 
several protocols. Moreover, one of our hypothesis was that sequence 
distribution linear symmetry plays a role on protease substract recoognition.
So we wanted to capture all sequences containing not only the entire interface
but neighbouring elements of sequence space.
<br>
*: The template can be any protein with sequence len > len_frags. It is just
a dummy molecule, so as to base itself upon and use Modellers optimization
algorithms and DOPE score function to produce "physically plausible" 
structures. These structures, upon docking initiallization, will probably
have all it's bond torsions randomized. But at, as with Deep Learnin, 
it all begins with a good method of initiallization. On a side note, I 
did implemented a more random way of doing this, with PeptideBuilder
library assembling fragments with phi-psi dihedrals randomized through
a [-pi,pi] uniform distribution.<br>
<br>
## Usage<br>
<br>
Using this script is simple. It expects the user too provide 3 required
inputs and one optional one. In orther to perform the entire routine 
(fragmentation and fragment modeling), it is mandatory to provide a 
multifasta file containing one or more polyprotein sequences (-p),
another multifasta for n interfaces (-i) and a receptor/template file (-r).
<br><br>
Like this:<br>
python3 run_example.py -p polyprot.fas -i interfaces.fas -r receptor_name \
-d docking_type (either VINA or AD4).<br>
<br>
It is noteworthy that this script will produce a ton of folders and directory
trees, soo keep in mind that the execution times will vary greatly upon:
sequence_length, num of sequences, number of different fragment lengths,
fragment lengths, type of preparation, and finally your computer specs.
<br><br>
## 1 minute tutorial:<br>
1- Create a (multi)fasta file containing one or more polyprotein aa sequence(s)<br>
2- Create another fasta for the interfaces (they all must have the same size,
so, in case of larger interfaces, croping is advised);<br>
3- Run the command below;<br>
4- Grab a cup of coffee and enjoy Modeller's verbose for several minutes.<br>
<br>
At the end all the models will be under a directory named after the original
fragment multifasta, produced by the Fragment Factory module.<br>
<br>
## ToDo:<br>
1- Implement a better logging mechanism.<br>
2- Reimplement the directory creation mechanism.<br>
3- Open to suggestions.<br>
