&run_type concat
#############################
#
# Inputs
#
############################
&input 
&name   2OH_on_Ti111-rlx.xyz
&type   xyz
&config 1
&end_input

&input 
&name   bromoisobutyrate_PBE_a2.92_Ecutoff400-evacc20.0-conf2-rlx.xyz
&type   xyz
&config 1
&end_input

&input 
&name   cp2k.restart
&type   cp2k_restart_file
&config 1
&end_input
############################
#
# Outputs
#
############################
&output
&name out.xyz
&type xyz
&end_output

#&output
#&name out.in
#&type cp2k_input
#&end_output


&machine
&name hpc
&end_machine
