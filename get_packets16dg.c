/* gcc -o get_packets16dg get_packets16dg.c -lm */
#include <stdio.h>
#include        <math.h>
/* for testing 16 antenna  configuration, 2015-10-27 */
#include "dpp_struct16.h"
#include "unp2.h"

void err_sys(const char* x);
void err_quit(const char* x);
void elapsed_time_ppp_(void);


/* the following are error_wrappers, used in the Unix Network Programming book */
void err_sys(const char* x)
{
  perror(x);
  exit(1);
}

void err_quit(const char* x)
{
  perror(x);
  exit(1);
}

/* A subprogram that tells me how long a process has been running */
void elapsed_time_ppp_(void)
{
  struct timeval otp;
  int rc;

  rc=gettimeofday(&otp, NULL);
  stime_.xsec = otp.tv_sec;
  stime_.xusec = otp.tv_usec;
  return;
}

/* Make these global, so that signal handler knows about them */
int                   sockfd;
struct sockaddr_in    servaddr, cliaddr;

void sig_handler(int signo)
{
  if (signo == SIGUSR1) 
  {
    printf("received SIGUSR1\n");
    close(sockfd);
    if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
      err_sys("socket error");
    bzero(&servaddr, sizeof(servaddr));
    servaddr.sin_family      = AF_INET;
    servaddr.sin_port        = htons(INPORT);
    if (inet_pton(AF_INET, argv[1], &servaddr.sin_addr) <= 0)
      err_quit("inet_pton error");
    if (bind(sockfd, (SA *) &servaddr, sizeof(servaddr)) < 0)
      err_quit("bind error");
  }
}

int
/* Changed so that the second argument is now the nth packet to be
   timed, jmm 2016-01-20  */
/* Changed so that there is no output, just time and accumulation
   testing, the second argument is back to the number of buffers */
main(int argc, char **argv)
{

  /*  static const char *IPADDR = "10.0.0.53";*/
  /*  static const char *IPADDR = "10.0.0.54";*/
  int                   sockfd;
  struct sockaddr_in    servaddr, cliaddr;
  int n, i, i1, j, ij, k, k1, k2, k3, slen=sizeof(cliaddr);
  int xsec0, xusec0;
  int nbuf, xcount=0, ycount=0, zcount=0; /*xcount is the count of packet dt's>5 msec */
  /* time test variables */
  double dtime, dtime2, dt1;
  double dtav, dtsig=0, dtmax=0, dtmin=10000.0;
  /* start denotes whether we've started recording, acnum is the accumulation number
     in a given second for the packet, gacc_num is the global accumulation number */
  unsigned int start, acc_num = 0, gacc_num = 0;
  unsigned int gacc_old = 0, acc_start = 0; /* the previous accumulation, and the start accumulation */
  unsigned int acc_rel = 0; /* This is the accumulation raltive to the start accum */
  unsigned int board_id = 0; /* The board_ID number */
  unsigned int pakno, xop; /* pakno is the packet number, xop is 1 for an X packet, 0 for a P */
  unsigned int pakno_sv; /* a saved packet number to check for non-contiguousness */

  unsigned char pak[PPAKBYTES];
  unsigned char hdr[HDRBYTES];
  unsigned char xpak[XPAKBYTES];
  int sxpak=sizeof(xpak);
  /* temporary buffer size */
  int bufsize_tmp=838860800, ibtmp, ibtmp1;
  struct timeval ctime;
  struct timeval tv; /*used for socket timeout */
  tv.tv_sec = 0;
  tv.tv_usec = 20000; /* 20 msec timeout */

  /* Temporary offset for X packet accumulation number */
  int offset_xxx=0;

  /* Set to max number of buffers (declared in dpp_struct.h) */
  if (argc>2)
    {
       nbuf = atoi(argv[2]);
    }

  time_t curtime;
  struct tm *loctime;
  char tbuffer[14];
  /* Get the current time.  */
  curtime = time (NULL);
  /* Convert it to local time representation.  */
  loctime = localtime (&curtime);
  strftime (tbuffer, 14, "%Y%m%d%H%M", loctime);
  printf("%s\n", tbuffer);

  char name1[31];

  strcpy(name1, "eovsa_packets16a_");
  strcat(name1, argv[1]);
  strcat(name1, ".txt");

  printf("argv[1] = %s\n", argv[1]);

  FILE *fl1 = fopen(name1, "w");
  if (!fl1)
    {
      fprintf(stderr, "Unable to open file %s", name1);
      return 0;
    }

  for (i=0; i<HDRBYTES; i++)
    {
      hdr[i]=0;
    }
  for (i=0; i<PPAKBYTES; i++)
    {
      pak[i]=0;
    }
  for (i=0; i<XPAKBYTES; i++)
    {
      xpak[i]=0;
    }
  for (i=0; i<NPAKCD; i++)
    {
      pakflag1_.flag[i]=255; /* a non-zero number means the packet has been processed */
    }
  for( i=0; i< NACCUMBUF; i++)
    {
      accumbuf.baccumno[i] = 0;
      accumbuf.npak[i] = 0;
    }

  signal(SIGUSR1, sig_handler)
  if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    err_sys("socket error");
  bzero(&servaddr, sizeof(servaddr));
  servaddr.sin_family      = AF_INET;
  servaddr.sin_port        = htons(INPORT);
  if (inet_pton(AF_INET, argv[1], &servaddr.sin_addr) <= 0)
    err_quit("inet_pton error");
  if (bind(sockfd, (SA *) &servaddr, sizeof(servaddr)) < 0)
    err_quit("bind error");
  //  if (setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, &bufsize_tmp, sizeof(bufsize_tmp)) < 0)
  //    err_quit("buffer size error");
  //  ibtmp = getsockopt(sockfd,SOL_SOCKET,SO_RCVBUF,&bufsize_tmp, &ibtmp1);
  //  printf("bufsize_tmp = %i\n", bufsize_tmp);
  //  setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(struct timeval)); /*set timeout option */

  k1 = 0; k2=0;
  dtime = 0; dtime2 = 0;
  start = 0;
  pakno = 0;
  for (i=0; i<nbuf; i++)
    {
      i1 = i % NBUF;
      for (j=0 ;j<NPACKET ;j++ )
        {
          elapsed_time_ppp_();
          xsec0 = stime_.xsec;
          xusec0 = stime_.xusec;
          //      usleep(2);
          //      bzero(&dppcorin_.correldata[k1], XPAKBYTES);
          n = recvfrom(sockfd, &dppcorin_.correldata[k1], XPAKBYTES, 0, (SA *) &cliaddr, &slen);
          if(n < 0)
            {
              printf("n = %i\n", n);
              err_quit("recvfrom error");
            }
          if(n == 0)
            {
              printf("n = %i\n", n);
              err_quit("recvfrom error");
            }
          /* See if I can extract the accumulation number from the buffer */
          /* acc_num = ntohs(*(uint16_t*)&dppcorin_.correldata[k1+13]); */
          elapsed_time_ppp_();
          dt1 = 1.0e6*(double) (stime_.xsec-xsec0) +
            (double) (stime_.xusec-xusec0);
          pakno_sv = pakno;
          memcpy(&gacc_num, &dppcorin_.correldata[ACCOFFSET+k1], 4);
          memcpy(&pakno, &dppcorin_.correldata[PAKNOOFFSET+k1], 2);
          memcpy(&board_id, &dppcorin_.correldata[BOARDIDOFFSET+k1], 2);
          memcpy(&acc_num, &dppcorin_.correldata[ACCINSECOFFSET+k1], 2);
          memcpy(&xop, &dppcorin_.correldata[XPOFFSET+k1], 2);
          acc_num = gacc_num % 50;
          ij = (j+i1*NPACKET);
          pkt_.acc_num[ij]=acc_num;
          pkt_.gacc_num[ij]=gacc_num;
          pkt_.pakno[ij]=pakno;
          pkt_.xop[ij]=xop;
          pkt_.board_id[ij]=board_id;
          pakflag1_.flag[ij]=0;
          if(start == 0)
            {
              acc_start = gacc_num;
              start = 1;
            }
          acc_rel = gacc_num-acc_start;
          if(acc_rel < NACCUMBUF)
            {
              accumbuf.baccumno[acc_rel] = gacc_num; /* accumulation number */
              accumbuf.npak[acc_rel]++; /* number of packets */
            }
          if(dt1 > 5000)
            {
              xcount++;
              //              printf("long count, %i %i %i %i\n", gacc_num, pakno, xop, board_id);
            }
          if(dt1 < 3)
            {
              ycount++;
            }
          if(dt1 > dtmax)
            {
              dtmax = dt1;
            }
          if(dt1 < dtmin)
            {
              dtmin = dt1;
            }
          if(dt1 == 0)
            {
              zcount++;
              //              printf("long count, %i %i %i %i\n", gacc_num, pakno, xop, board_id);
            }
          if((gacc_num - gacc_old) > 2)
            {
              printf("gacc_num = %i , gacc_old = %i\n", gacc_num, gacc_old);
            }
          gacc_old = gacc_num;
          dtime = dtime+dt1;
          dtime2 = dtime2+dt1*dt1;
          k1 = (k1+XPAKBYTES) % NCD;
          k2++; /* k2 will be the number of packets processed */
        }
    }

  /* how long did it take */
  dtav = dtime/(k2);
  dtsig = sqrt((dtime2-k2*dtav*dtav)/(k2-1));

  printf("Time to process each packet in microseconds: %f\n", dtav);
  printf("Standard deviation of time to process each packet in microseconds: %f\n", dtsig);
  printf("Max time to process each packet in microseconds: %f\n", dtmax);
  printf("Min time to process each packet in microseconds: %f\n", dtmin);
  printf("Number of packets that took longer than 5 msec: %i\n", xcount);
  printf("Number of packets that took less than 3 musec: %i\n", ycount);
  printf("Number of packets that took 0 musec: %i\n", zcount);
  printf("Total number of packets %i\n", k2);
  printf("Start accumulation: %i\n", acc_start);
  printf("End accumulation: %i\n", acc_rel+acc_start);
  printf("Accumulation Number, Npackets:\n");
  fprintf(fl1, "Start accumulation: %i\n", acc_start);
  fprintf(fl1, "End accumulation: %i\n", acc_rel+acc_start);
  fprintf(fl1, "Accumulation Number, Npackets:\n");
  for (j=0; j<NACCUMBUF ;j++)
    {
     if(accumbuf.baccumno[j] > 0 && accumbuf.npak[j] < 3072)
       //       if(accumbuf.baccumno[j] > 0)
        {
          printf("%i  %i\n", accumbuf.baccumno[j], accumbuf.npak[j]);
          fprintf(fl1, "%i  %i\n", accumbuf.baccumno[j], accumbuf.npak[j]);
        }
    }

  fprintf(fl1, "Time to process each packet in microseconds: %f\n", dtav);
  fprintf(fl1, "Standard deviation of time to process each packet in microseconds: %f\n", dtsig);
  fprintf(fl1, "Max time to process each packet in microseconds: %f\n", dtmax);
  fprintf(fl1, "Min time to process each packet in microseconds: %f\n", dtmin);
  fprintf(fl1, "Number of packets that took longer than 5 msec: %i\n", xcount);
  fprintf(fl1, "Number of packets that took less than 3 musec: %i\n", ycount);
  fprintf(fl1, "Number of packets that took 0 musec: %i\n", zcount);
  fprintf(fl1, "Total number of packets %i\n", k2);
  fclose(fl1);
  close(sockfd);
  return 0;
}
