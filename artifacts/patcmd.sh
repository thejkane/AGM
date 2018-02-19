#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+25507-212s.ap2 
#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+30436-212s.ap2
#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+30088-212s.ap2

#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+34180-212s.ap2
#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+71163-212s.ap2
#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times /home/users/p02119/development/pbgl2/build/bin/mis_family+87956-212s.ap2
#pat_report /home/users/p02119/development/pbgl2/build/bin/mis_family+87956-212s.ap2

#pat_report -s th=ALL -O affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times $1
pat_report -s th=ALL -O hwpc,affinity,load_imbalance_thread,profile_th_pe,profile_pe_th,thread_times $1

#-O hwpc

