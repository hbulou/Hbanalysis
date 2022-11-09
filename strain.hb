&run_type strain
&system strain

&input 
&name   essai-1.restart
&type   cp2k_restart_file
&end_input

&input 
&name   essai-pos-1.xyz
&type   xyz
&end_input

&transformation
  &strain TRUE
  &eps
  0.0 0.0 0.0
  0.0 0.0 0.0
  0.0 0.0 0.01
  &end_strain
&end_transformation


&output 
&name   test.in
&type   cp2k_restart_file
&end_output

