#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "heidi.h"
#include "dynamics.h"
#include "io.h"


Environment env;


void help()
{
    fprintf(stderr, "heidi <filename>\n");
    exit(2);
}

int main(int argc, char **argv)
{
    int i;

    if (argc < 2) help();

    static struct option long_options[] = {
        {"state", required_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    init_env(&env);

    /*========================================================================
     * Process the command line flags
     *======================================================================*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "s:vlo:",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                break;

            case 's': 
                strncpy(env.io.state_file, optarg, sizeof(env.io.state_file)-1);
                env.io.state_file[sizeof(env.io.state_file)-1] = '\0';
                break;
            case 'v': env.io.verbosity++;                      break;

            case 'h': help(); break;
            case '?': break;
        }
    }

    env.t_now  = 0;
    env.t_end  = SECONDS(1e-6);
    env.t_step = 0;

    VL(2)
    {
        fprintf(env.io.err, "%e %e\n", env.t_end, SECONDS(1));
#define print_sizeof(type) do { fprintf(env.io.err, "%-30s = %4ld bytes\n", "sizeof(" #type ")", sizeof(type)); } while(0)
        print_sizeof(Particle);
        print_sizeof(ParticleTypeData);
        print_sizeof(Environment);
        print_sizeof(Constants);
    }

    load_particles(argv[optind]);

    for all_particles(i) P_dt(i) = SECONDS(2e-10);


    //while (env.t_now < env.t_end)
    while (P_xx(0) < 5e35)
    {
        write_state();
        Time dt = step_particles();

        env.t_now += dt;
        env.t_step++;

        if ((env.t_step % 100) == 0)
            VL(2) fprintf(env.io.err, ".");
    }


    io_close(env.io.state_fp);
    io_close(env.io.stats_fp);

    return 0;
}

