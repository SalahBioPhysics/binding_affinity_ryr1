import sys
import os
import time
import os.path
from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call # This is needed to submit jobs

"""
This code is used to submit mcce jobs since we have 50 systems 
"""
 
pdb_files = '/home/guest/ryr1_energyCalculations/input_data/activation_core/new_model/'
destination_runs = '/home/guest/ryr1_energyCalculations/calculations/atp_large_vdw/'


def change_runprm(runprm,prot,str1,str2,str3,str4):
    lines_ = open(runprm).readlines()
    for line in lines_:
        a = line
        line = line.split()
        #print line[-1]
        create_new = runprm+'2'
        with open(create_new,'a') as prm:
            if len(line) != 0:
                if line[-1] == '(INPDB)':
                    newLine = prot+'   '+line[-1]+'\n'
                    prm.write(newLine)
                elif line[-1] == '(DO_PREMCCE)':
                    newLine = str1+'    '+line[-1]+'\n'
                    prm.write(newLine)
                elif line[-1] == '(DO_ROTAMERS)':
                    newLine = str2+'    '+line[-1]+'\n'
                    prm.write(newLine)
                elif line[-1] == '(DO_ENERGY)':
                    newLine = str3+'    '+line[-1]+'\n'
                    prm.write(newLine)
                elif line[-1] == '(DO_MONTE)':
                    newLine = str4+'    '+line[-1]+'\n'
                    prm.write(newLine)
                else:
                    prm.write(a)

            prm.close()
    sys_call = 'mv ' + create_new + ' ' + runprm
    os.system(sys_call) 
 

def head3(head3_file,mydirectory_2):
    template = '{0:6}{1:15}{2:2}{3:6}{4:10}{5:3}{6:6}{7:3}{8:3}{9:8}{10:9}{11:6}{12:9}{13:9}{14:6}{15:11}\n'
    lines_ = open(head3_file).readlines()
    create_new = mydirectory_2+'head3.lst'
    with open(create_new,'w') as fp:
        for line in lines_:
            fields = line.split()
            if fields[1] == 'ATP-4A0001_001':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],fields[2],float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
                #fp.write(template.format(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],'21.960',fields[11],fields[12],fields[13],fields[14],fields[15]))
            elif fields[1] == 'ATP-AA0001_002':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-AA0001_003':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-AA0001_004':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-BA0001_005':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-BA0001_006':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-BA0001_007':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            elif fields[1] == 'ATP-CA0001_008':
                fp.write("%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n" % (int(fields[0]),fields[1],'t',float(fields[3]),float(fields[4]),float(fields[5]),float(fields[6]),int(fields[7]),int(fields[8]),00.000,00.000,float(fields[11]),float(fields[12]),float(fields[13]),float(fields[14]),fields[15]))
            else:
                    fp.write(line)
        fp.close()





for i in range(1,51):
    #onlyfiles = [f for f in os.listdir(pdb_files+topDir) if os.path.isfile(os.path.join(pdb_files+topDir, f))]
    #pdbfile = topDir.split('-')
    mydirectory = destination_runs + "frame_" +str(i).zfill(2)+'/'
    if not os.path.exists(mydirectory):
        sys_call = 'mkdir ' + mydirectory
        os.system(sys_call)

    #mydirectory_2 = destination_runs + "frame_" +str(i).zfill(2)+'/ATP01_new/'
    #if not os.path.exists(mydirectory_2):
        #sys_call = 'mkdir ' + mydirectory_2
        #os.system(sys_call)
    


    # move submit.sh and run.prm
    sys_call = 'cp ' + destination_runs+'submit.sh '+mydirectory
    os.system(sys_call)
    
    sys_call = 'cp ' + destination_runs+'run.prm '+mydirectory
    os.system(sys_call)
    # move the pdb files
    sys_call = 'cp ' + pdb_files + str(i).zfill(2)+'.pdb ' + mydirectory
    os.system(sys_call)

    # move energies.opp
    #sys_call = 'cp ' + mydirectory+'energies.opp '+mydirectory_2
    #os.system(sys_call)

    # move energies.opp
    #sys_call = 'cp ' + mydirectory+'head3.lst '+mydirectory_2
    #os.system(sys_call)

    
    #head3_file = mydirectory+'/head3.lst'
    #head3(head3_file,mydirectory_2)


    
    
    # move head3.lst, step2out, and energies.opp 
    #sys_call = 'cp ' + destination_runs_head3 + "frame_" +str(i).zfill(2)+'/energies.opp '+mydirectory
    #os.system(sys_call)

    # Change the (INPUT) param in run.prm to say the correct file name (str(i).zfill(2)+'.pdb ')
    # 1. runprm: the run.prm file that we will edit
    # 2. prot: the name and location of the PDB file
    # 3. str1, str2, str3 and str4 are true (T) or false (F) flags 
    

    runprm = mydirectory+'run.prm'
    prot   = mydirectory+str(i).zfill(2)+'.pdb'
    change_runprm(runprm,prot,"t","t","t","f")

    # Send to submit 
    os.chdir(mydirectory)
    qsub_call = "qsub %s"
    print ('doing ' + mydirectory)
    call(qsub_call % "submit.sh", shell=True)
    # Wait for 5 seconds, so we don't overwork the queue system
    time.sleep(2)
 
print ('Done')

