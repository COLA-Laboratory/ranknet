//
// Created by gylai on 2021/8/25.
//

#include "weights_by_har.h"

int main(int argc, char* argv[])
{
    char input_file[MAX_CHAR_SIZE];
    char har_file[MAX_CHAR_SIZE];
    char sh_file[MAX_CHAR_SIZE];
    char *para_name;

    int nsize, ndimension;

    for (int i = 1; i+1 < argc; ++i)
    {
        para_name = argv[i]+1;
        if (!strcmp(para_name, "d"))
        {
            ndimension = atoi(argv[i+1]);
        }
        else if (!strcmp(para_name, "n"))
        {
            nsize = atoi(argv[i+1]);
        }
    }

    double boundaries[ndimension][2];


    sprintf(input_file, "../../hitandrun/input_%dd_%d.txt", ndimension, nsize);
    sprintf(sh_file, "../../hitandrun/run_%dd_%d.sh", ndimension, nsize);
    sprintf(har_file, "../../UniformWeights/%dd_%d.txt", ndimension, nsize);

    //initialize the boundaries
    for (int i = 0; i < ndimension; ++i)
    {
        boundaries[i][0] = 0.0;
        boundaries[i][1] = 1.0;
    }
    std::ofstream out;
    out.open(input_file);
    for (int i = 0; i < ndimension; ++i)
    {
        for (int j = 0; j < ndimension; ++j)
        {
            out<<(i==j)<<" ";
        }
        out<<">= "<<boundaries[i][0]<<std::endl;
        for (int j = 0; j < ndimension; ++j)
        {
            out<<(i==j)<<" ";
        }
        out<<"<= "<<boundaries[i][1]<<std::endl;
    }
    for (int i = 0; i < ndimension; ++i)
    {
        out<<1<<" ";
    }
    out<<"= "<<1;
    out.flush();
    out.close();
    
    out.open(sh_file);
    char shell_command[MAX_CHAR_SIZE];
    sprintf(shell_command, "java -jar ../../polyrun-1.0.0-jar-with-dependencies.jar -i %s -n %d > %s\n"
                   "wait\n"
                   "exit",
            input_file, nsize, har_file);
    out<<shell_command;
    out.flush();
    out.close();

    char giveRight[MAX_CHAR_SIZE];
    sprintf(giveRight, "chmod 777 %s", sh_file);
    system(giveRight);

    pid_t status = system(sh_file);
    if (-1 == status)
    {
        std::cerr<<"<error> system error!\n";
    }
    else
    {
        printf("<info> generate %d weights in %d dimension(s) via hitandrun.\nsee in %s\n",
               nsize, ndimension, har_file);
    }

    return 0;
}




