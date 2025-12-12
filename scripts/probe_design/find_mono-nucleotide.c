/* Example usage:
 * gcc find_mono-nucleotide.c
 * cat /media/Scratch_SSD_Voyager/picheng/GenomeDB/mm39/raw_fasta/mm39.fa | ./a.out >mm39_MS_location.mono-nucleotide_10bp-max.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define SEQNAME_MAX BUFSIZ
#define REPORT if(len_repeat>9) printf("010\t%s\t%d\t%d\t%dx%c\n",name,pos-len_repeat,pos,len_repeat,prev_c); len_repeat=0;
int main(int argc,char** argv)
    {
    int c;
    int pos=0;
    int prev_c=-1;
    char name[SEQNAME_MAX];
    name[0]=0;
	int len_repeat=0;
    for(;;)
            {
            switch((c=fgetc(stdin)))
                {
				case EOF: return EXIT_SUCCESS;
				case '>':
				    {
				    int space=0;
				    int name_length=0;
				    REPORT;
				    name[0]=0;
				    pos=0;
				    while((c=fgetc(stdin))!=EOF && c!='\n')
				        {    
				        if(space) continue;
				        if(isspace(c)) { space=1; continue;}
				        name[name_length++]=c;
				        }
				    name[name_length]=0;
				    prev_c=-1;
				    len_repeat=0;
				    break;
				    }
				case '\n':case '\r':case ' ':break;
				case 'a': case 'A':
				case 't': case 'T':
				case 'g': case 'G':
				case 'c': case 'C':
					{
					c= toupper(c);
				    if(prev_c==c)
				        {
				       	++len_repeat;
				        }
				    else
				    	{
				    	REPORT;
				    	}
				    prev_c=c;
				    ++pos;
				    break;
				    }
                default:prev_c=c; ++pos;break;
                }
        }
    return EXIT_SUCCESS;
    }
