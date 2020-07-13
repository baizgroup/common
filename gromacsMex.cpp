/* this is a modified gromacs xtc-parser in combination with mex/matlab
 * Created 2012-2014
 * Last modified: 4.9.2013
 * Version 1.0
 */

/*
 * This code is partly based on the XTC loader module of VMD (http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/Gromacs_8h-source.html)
 * and the GROMACS libxdrf (ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.5.tar.gz)
 *
 */

#include <matrix.h> // from tutorial

#include <mex.h>   
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
// #include "endianswap.h" // not required?


// Error codes for mdio_errno
#define MDIO_SUCCESS		0
#define MDIO_BADFORMAT		1
#define MDIO_EOF			2
#define MDIO_BADPARAMS		3
#define MDIO_IOERROR		4
#define MDIO_BADPRECISION	5
#define MDIO_BADMALLOC		6
#define MDIO_CANTOPEN		7
#define MDIO_BADEXTENSION	8
#define MDIO_UNKNOWNFMT		9
#define MDIO_CANTCLOSE		10
#define MDIO_WRONGFORMAT	11
#define MDIO_SIZEERROR		12
#define MDIO_UNKNOWNERROR	1000

#define MDIO_READ	0
#define MDIO_WRITE	1

#define MDIO_MAX_ERRVAL		11

static int mdio_errcode;	// Last error code

#define TRX_MAGIC	1993	// Magic number for .trX files
#define XTC_MAGIC	1995	// Magic number for .xtc files
#define MAX_GRO_LINE	500	// Maximum line length of .gro files
#define MAX_G96_LINE	500	// Maximum line length of .g96 files
#define MAX_TRX_TITLE	80	// Maximum length of a title in .trX
#define MAX_MDIO_TITLE	80	// Maximum supported title length
#define ANGS_PER_NM		10	// Unit conversion factor
#define ANGS2_PER_NM2	100	// Unit conversion factor


// All the supported file types and their respective extensions
#define MDFMT_GRO		1
#define MDFMT_TRR		2
#define MDFMT_G96		3
#define MDFMT_TRJ		4
#define MDFMT_XTC		5

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(xtc_magicints) / sizeof(*xtc_magicints))

typedef struct {
	int version;		// File version number
	char title[MAX_TRX_TITLE + 1];	// File title
	int ir_size;
	int e_size;
	int box_size;
	int vir_size;
	int pres_size;
	int top_size;
	int sym_size;
	int x_size;		// Positions of atoms
	int v_size;		// Velocities of atoms
	int f_size;
	int natoms;		// Number of atoms in the system
	int step;
	int nre;
	float t;
	float lambda;
} trx_hdr;


// A generic i/o structure that contains information about the
// file itself and the input/output state
typedef struct {
	FILE *	f;	// Pointer to the file
	int	fmt;	// The file format
	int	prec;	// Real number precision
	int	rev;	// Reverse endiannism?
	trx_hdr * trx;	// Trx files require a great deal more
			// header data to be stored.
} md_file;

// A format-independent structure to hold unit cell data
typedef struct {
  float A, B, C, alpha, beta, gamma;
} md_box;



// Timestep information
typedef struct {
	float *pos;	// Position array (3 * natoms)
	//float *vel;	// Velocity array ** (VMD doesn't use this) **
	//float *f;	// Force array ** (VMD doesn't use this) **
	//float *box;	// Computational box ** (VMD doesn't use this) **
	int natoms;	// Number of atoms
	int step;	// Simulation step
	float time;	// Time of simulation
  	md_box *box;
} md_ts;


/*
 vars from tutorial
*/
/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif
/*
 end vars from tutorial
*/

/***********************************************************************

FUNCTIONS 

***********************************************************************/

static int mdio_seterror(int);

// integer table used in decompression
static int xtc_magicints[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0,8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
	80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
	1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384,
	20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031, 131072,
	165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255,
	1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304,
	5284491, 6658042, 8388607, 10568983, 13316085, 16777216 };


// reads bits from a buffer.    
static int xtc_receivebits(int *buf, int nbits) {

	//printf("#IN %i #",  buf[2]);
	int cnt, num; 
	unsigned int lastbits, lastbyte;
	unsigned char * cbuf;
	int mask = (1 << nbits) -1;
	// sizeof buf -> 4, as 4 elements of an array
	cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);  // equivalent: unsigned char* cbuf = (unsigned char*) &buf[3]; 
	cnt = buf[0];
	lastbits = (unsigned int) buf[1];
	lastbyte = (unsigned int) buf[2];
	num = 0;
	
	while (nbits >= 8) {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (nbits - 8);
		nbits -=8;
	}
		
	if (nbits > 0) {
		if (lastbits < (unsigned int)nbits) {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
			
		}
		lastbits -= nbits;
		num |= (lastbyte >> lastbits) & ((1 << nbits) -1);
	}

	num &= mask;
	buf[0] = cnt;
	buf[1] = lastbits;
	buf[2] = lastbyte;
	return num; 
}


// decompresses small integers from the buffer
static void xtc_receiveints(int *buf, const int nints, int nbits, unsigned int *sizes, int *nums) {
	int bytes[32];
	int i, j, nbytes,p,  num;

	bytes[1] = bytes[2] = bytes[3] = 0;
	nbytes = 0;
	i=0;	
	do{	
		bytes[i] =0;
		i++;
	}while(i<32);

	while (nbits > 8) {
		bytes[nbytes++] = xtc_receivebits(buf, 8);
		nbits -= 8;
	}
	
	if (nbits > 0) {
		bytes[nbytes++] = xtc_receivebits(buf, nbits);
	}

	for (i = nints-1; i > 0; i--) {
		num = 0;

		for (j = nbytes-1; j >=0; j--) {
		
		
			num = (num << 8) | bytes[j];
			p = num /sizes[i]; // bug from chromacs? num:sizes,  -341664787 : 1358 = 2911121 // wrong result -> bug?? in matlab i get the correct value
			// reason -> sizes is an array and change the division result..if i copy sizes to int tmp variable and then i divide i get the correct result
			// printf("%i: %i, %i \n", (j+1), num, sizes[i]);
			bytes[j] = p;
			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

// xtc_data() - reads a specific amount of data from an xtc
// file using the xdr format.
static int xtc_data(md_file *mf, char *buf, int len) {
	size_t slen = (size_t)len;
	if (!mf || len < 1) return mdio_seterror(MDIO_BADPARAMS);
	
	if (buf) {
		if (fread(buf, 1, slen, mf->f) != slen) {
			if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
			if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
			else return mdio_seterror(MDIO_UNKNOWNERROR);
		}
		if (len % 4) {
			if (fseek(mf->f, 4 - (len % 4), SEEK_CUR)) {
				if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
				if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
				else return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
	}else{
		int newlen;
		newlen = len;
		if (len % 4) newlen += (4 - (len % 4));
		if (fseek(mf->f, newlen, SEEK_CUR)) {
			if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
			if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
			else return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}
	return len;
}

// returns the number of bits in the binary expansion of
// the given integer.
static int xtc_sizeofint(int size) {
	unsigned int num = 1;
  unsigned int ssize = (unsigned int)size;
	int nbits = 0;

	while (ssize >= num && nbits < 32) {
		nbits++;
		num <<= 1;
	}
	return nbits;
}

// calculates the number of bits a set of integers, when compressed,
// will take up.
static int xtc_sizeofints(int nints, unsigned int *sizes) {
	int i;
  unsigned int num;
	unsigned int nbytes, nbits, bytes[32], bytecnt, tmp;
	nbytes = 1;
	bytes[0] = 1;
	nbits = 0;
	for (i=0; i < nints; i++) {	
		tmp = 0;
		for (bytecnt = 0; bytecnt < nbytes; bytecnt++) {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}

		while (tmp != 0) {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;

		}
		nbytes = bytecnt;
	}
	num = 1;
	nbytes--;

	while (bytes[nbytes] >= num) {
		nbits++;
		num *= 2;
		
	}

	return nbits + nbytes * 8;
}


// Sets the error code and returns an appropriate return value
// for the calling function to return to its parent
static int mdio_seterror(int code) {
	mdio_errcode = code;
	return code ? -1 : 0;
}

// Converts box basis vectors to A, B, C, alpha, beta, and gamma.  
// Stores values in md_box struct, which should be allocated before calling
// this function.
static int mdio_readbox(md_box *box, float *x, float *y, float *z) {
  float A, B, C;

  if (!box) {
    return mdio_seterror(MDIO_BADPARAMS);
  }

  // A, B, C are the lengths of the x, y, z vectors, respectively

  A = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) * ANGS_PER_NM;
  B = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] ) * ANGS_PER_NM;
  C = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] ) * ANGS_PER_NM;

  if ((A<=0) || (B<=0) || (C<=0)) {
    /* Use zero-length box size and set angles to 90. */
    box->A = box->B = box->C = 0;
    box->alpha = box->beta = box->gamma = 90;
  } else {
    box->A = A;
    box->B = B;
    box->C = C;
    // gamma, beta, alpha are the angles between the x & y, x & z, y & z
    // vectors, respectively
    box->gamma = acos( (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])*ANGS2_PER_NM2/(A*B) ) * 90.0/M_PI_2;
    box->beta = acos( (x[0]*z[0]+x[1]*z[1]+x[2]*z[2])*ANGS2_PER_NM2/(A*C) ) * 90.0/M_PI_2;
    box->alpha = acos( (y[0]*z[0]+y[1]*z[1]+y[2]*z[2])*ANGS2_PER_NM2/(B*C) ) * 90.0/M_PI_2; 
  }

  return mdio_seterror(MDIO_SUCCESS);
}

// xtc_int() - reads an integer from an xtc file
static int xtc_int(md_file *mf, int *i) {
	unsigned char c[4];

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);
    
	// sanity check.
    if (sizeof(int) != 4) return mdio_seterror(MDIO_SIZEERROR);

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
		else if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
		else return mdio_seterror(MDIO_UNKNOWNERROR);
	}
	
	if (i) *i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
	return mdio_seterror(MDIO_SUCCESS);
}


// xtc_float() - reads a float from an xtc file
static int xtc_float(md_file *mf, float *f) {
	unsigned char c[4];
	int i;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
		else if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
		else return mdio_seterror(MDIO_UNKNOWNERROR);
	}

	if (f) {
		// By reading the number in as an integer and then
		// copying it to a floating point number we can
		// ensure proper endianness
		i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
		memcpy(f, &i, 4); // copy data to md_file structure md_ts (i.e. ts->pos)
	}
	return mdio_seterror(MDIO_SUCCESS);
}

// function that actually reads and writes compressed coordinates    
static int xtc_3dfcoord(md_file *mf, float *fp, int *size, float *precision, int frame) {


	/* original code
	static int *ip = NULL;
	static int oldsize;
	static int *buf;
*/
	// my code, deleted "static" to make sure that the program can be called more than once without interruptions
	 int *ip = NULL;
	 int oldsize;
	 int *buf;
	// end my code

	int minint[3], maxint[3], *lip;
	int smallidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int flag, k;
	int small, smaller, i, is_smaller, run;
	float *lfp;
	int tmp, *thiscoord,  prevcoord[3];

	int bufsize, lsize;
	unsigned int bitsize;
	float inv_precision;

	if (xtc_int(mf, &lsize) < 0) return -1;
	if (*size != 0 && lsize != *size) return mdio_seterror(MDIO_BADFORMAT);

	*size = lsize;
	size3 = *size * 3;

	if (*size <= 9) {
		for (i = 0; i < *size; i++) {
			if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
		}
		return *size;
	}

	xtc_float(mf, precision);

	if (ip == NULL) {
		ip = (int *)mxMalloc(size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)mxMalloc(bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	} else if (*size > oldsize) {
		ip = (int *)realloc(ip, size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)realloc(buf, bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	}

	buf[0] = buf[1] = buf[2] = 0;
	xtc_int(mf, &(minint[0]));
	xtc_int(mf, &(minint[1]));
	xtc_int(mf, &(minint[2]));
	xtc_int(mf, &(maxint[0]));
	xtc_int(mf, &(maxint[1]));
	xtc_int(mf, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
		bitsizeint[0] = xtc_sizeofint(sizeint[0]);
		bitsizeint[1] = xtc_sizeofint(sizeint[1]);
		bitsizeint[2] = xtc_sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	} else {
		bitsize = xtc_sizeofints(3, sizeint);
	}
	xtc_int(mf, &smallidx);
	smaller = xtc_magicints[FIRSTIDX > smallidx - 1 ? FIRSTIDX : smallidx - 1] / 2;
	small = xtc_magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;

	/* buf[0] holds the length in bytes */
	
	if (xtc_int(mf, &(buf[0])) < 0) return -1;
	
	/* result: return -1, y? */
	if (xtc_data(mf, (char *) &buf[3], (int) buf[0]) < 0) return -1;
	
	buf[0] = buf[1] = buf[2] = 0;
	lfp = fp;
	inv_precision = 1.0f / (*precision);
	run = 0;
	i = 0;
	lip = ip;

	while (i < lsize) {
		thiscoord = (int *)(lip) + i * 3; // HD: effects that thiscoord starts always with no value in the loop??
		if (bitsize == 0) {
			// hd: in this case this code will be never loaded
			thiscoord[0] = xtc_receivebits(buf, bitsizeint[0]);
			thiscoord[1] = xtc_receivebits(buf, bitsizeint[1]);
			thiscoord[2] = xtc_receivebits(buf, bitsizeint[2]);		
		} else {	
			xtc_receiveints(buf, 3, bitsize, sizeint, thiscoord);
		}

		i++;

		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];
/*

char c; 
scanf("%c",&c); 
while(getchar() != '\n');
*/
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];

		flag = xtc_receivebits(buf, 1);
		is_smaller = 0;

		if (flag == 1) {
			run = xtc_receivebits(buf, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
	
		if (run > 0) {
			thiscoord += 3; // HD note: just effects that all elements in the array thiscoord gets the value 0
			for (k = 0; k < run; k+=3) {
				
				xtc_receiveints(buf, 3, smallidx, sizesmall, thiscoord);
				i++;
				thiscoord[0] += prevcoord[0] - small;
				thiscoord[1] += prevcoord[1] - small;
				thiscoord[2] += prevcoord[2] - small;
			
				
				if (k == 0) {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
				
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					
					*lfp++ = prevcoord[0] * inv_precision;
					*lfp++ = prevcoord[1] * inv_precision;
					*lfp++ = prevcoord[2] * inv_precision;
					
				}else{
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
					

				}
			
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;
			} // loop for
				
		} else {
			*lfp++ = thiscoord[0] * inv_precision;
			*lfp++ = thiscoord[1] * inv_precision;
			*lfp++ = thiscoord[2] * inv_precision;	
		}

		smallidx += is_smaller;
		// printf("%i, ", smallidx);
		if (is_smaller < 0) {
			small = smaller;
			if (smallidx > FIRSTIDX) {
				smaller = xtc_magicints[smallidx - 1] /2;
			} else {
				smaller = 0;
			}
		} else if (is_smaller > 0) {
			smaller = small;
			small = xtc_magicints[smallidx] / 2;
		}
		
		sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;

	} // loop

	return 1;
}


//static int xtc_timestep(md_file *mf, md_ts *ts) { // orignal version
static int xtc_timestep(md_file *mf, md_ts *ts, int frame, float data[], int counter, int numberAtoms, 
        int step[], float time[], float retPrecision[], float box[], int startFrame, int showmaxframes ) {

	float f, x[3], y[3], z[3];
	float precision=0;
	int size = 0;
	int n;
    int k=0;
	int j=0;
    int row=0;
    float floatTmp;
    int intTmp;
    

	// Check magic number
	if (xtc_int(mf, &n) < 0) return 1995;
    if (n != XTC_MAGIC) return 1995;

	// Get number of atoms
	if (xtc_int(mf, &n) < 0) return -1;
	ts->natoms = n;
   
	// Get the simulation step
	if (xtc_int(mf, &n) < 0) return -1;
	ts->step = n;

	// Get the time value
	if (xtc_float(mf, &f) < 0) return -1;
	ts->time = f;

	// Read the basis vectors of the box
  	if ( (xtc_float(mf, &x[0]) < 0) ||
      	 (xtc_float(mf, &x[1]) < 0) ||
      	 (xtc_float(mf, &x[2]) < 0) ||
      	 (xtc_float(mf, &y[0]) < 0) ||
      	 (xtc_float(mf, &y[1]) < 0) ||
      	 (xtc_float(mf, &y[2]) < 0) ||
      	 (xtc_float(mf, &z[0]) < 0) ||
      	 (xtc_float(mf, &z[1]) < 0) ||
      	 (xtc_float(mf, &z[2]) < 0) )
    	return -1;

  // Allocate the box and convert the vectors.
  ts->box = (md_box *) mxMalloc(sizeof(md_box));

  if (mdio_readbox(ts->box, x, y, z) < 0) {
    free(ts->box);
    ts->box = NULL;
    return -1;
  }
	
	ts->pos = (float *) mxMalloc(sizeof(float) * 3 * ts->natoms);
	if (!ts->pos) return mdio_seterror(MDIO_BADMALLOC);
	n = xtc_3dfcoord(mf, ts->pos, &size, &precision, frame);

        
	/* scaling... */
	for (n = 0; n < ts->natoms * 3; n++)
		ts->pos[n] *= ANGS_PER_NM;
    
    
    row=1;
    counter=0;
    step[frame] = ts->step;
    time[frame] = ts->time;
    retPrecision[frame] = precision;
    
    box[frame] = ts->box->A/ANGS_PER_NM;
    box[frame+showmaxframes] = ts->box->B/ANGS_PER_NM;
    box[frame+showmaxframes*2] = ts->box->C/ANGS_PER_NM;

    
    do {
      // uncomment
       // if(ts->pos[j] < 0) {
       //     break;
       // }

       if(row > 3){
              row=1;
              counter++;
       }

       floatTmp = ts->pos[j]/ANGS_PER_NM;
       if(row==1){
            data[counter+(frame*numberAtoms*3)] = floatTmp;
       }else if(row==2){
            data[counter+numberAtoms+(frame*numberAtoms*3)] = floatTmp;
       }else if(row==3){
            data[counter+numberAtoms*2+(frame*numberAtoms*3)] = floatTmp;
       }
       
       row++;
       j++;
    } while (j < numberAtoms*3);  

         
    return mdio_seterror(MDIO_SUCCESS);
}

 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int chkMagicN=0;
    int chknAtoms=0;
    int countFrame=1;
    long tellPos=0;
    int endOfFile=0;
    int numberAtoms=0;
    FILE * pFile;
    
    int showmaxframes=1;
	int frame=0;
	int breakVal=0;
    int counter =0;

    float *data2;
    int *dataInt;
    int *step;
    float *time;
    float *precision;
    float *box;
    
    int startFrame=0;
    int endFrame=0;

    mxArray *xData;
    char *filename;
    
    mwSize dim[3] = {0,3,0}; // first dimension zero will be replaced by argument 
    mwSize dim2[1] = {0};    // steps, time, precision
    mwSize dim3[2] = {0, 3}; // box
    
 	md_file *mf;
 	md_ts *ts;
	mf = (md_file *) mxMalloc(sizeof(md_file));
	ts = (md_ts *) mxMalloc(sizeof(md_ts));
   
    /*
     * get parameters
     */
    xData = (mxArray *)prhs[0];
    startFrame = mxGetScalar(xData);
    
    xData = (mxArray *)prhs[1];
    endFrame = mxGetScalar(xData);

    filename = mxArrayToString(prhs[2]);      
    mf->f = fopen(filename, "rb"); 
    
    // if file doesnt exist
    if(mf->f == false){
        printf("\nFile doesn't exist! \n");
        return;
    }
   
     // check if file is a xtc-file
     xtc_int(mf, &chkMagicN); 
     if(chkMagicN != 1995) {
        printf("\nFile is not a XTC-File! \n");
        return;     
     }

     // get numberOfAtoms for double chk on data records
     xtc_int(mf, &numberAtoms); // number of atoms
       
    // rewind
     rewind (mf->f);
    
    //
    // count number of all frames
    //
    while (!feof(mf->f)) {
        xtc_int(mf, &chkMagicN); 
        if(chkMagicN==1995) {
            xtc_int(mf, &chknAtoms);
            if(chknAtoms == numberAtoms) {              
                if(countFrame == startFrame) {    
                    tellPos = ftell(mf->f)-8;  // remember pos
                }
                
                countFrame++;
                if(endFrame != 0 &&  (countFrame-1) > endFrame) {
                    break;
                }               
                fseek ( mf->f , -4 , SEEK_CUR );
            }
        }
        
        if(feof(mf->f))
            endOfFile = 1;
        
    }
    
    // chk arguments
    if(startFrame > (countFrame-1) && endOfFile) {
        printf("Startframe does not exist! %i frames are available. \n", countFrame-1);
        return;
    }

    if(endFrame == 0 || endOfFile) {
        endFrame = countFrame-1;   
    }

    fseek ( mf->f , tellPos , SEEK_SET ); // move pointer back cause of 1995
    showmaxframes = endFrame-startFrame+1; 
    
    dim[0] = numberAtoms; 
    dim[2] = showmaxframes;
    dim2[0] = showmaxframes;
    dim3[0] = showmaxframes;
    
    // allocate memory
    data2 = (float *)mxMalloc(numberAtoms*3*showmaxframes*sizeof(*data2)); 
    dataInt = (int *)mxMalloc(showmaxframes*sizeof(*dataInt));
    step =  (int *)mxMalloc(showmaxframes*sizeof(*step));
    time =  (float *)mxMalloc(showmaxframes*sizeof(*time));
    precision =  (float *)mxMalloc(showmaxframes*sizeof(*precision));
    box =  (float *) mxMalloc(3*showmaxframes*sizeof(*box));
 
    /*
     * Read file
     */
    
   if(mf->f != NULL) {
     // coords
     plhs[0] = mxCreateNumericArray(3, dim, mxSINGLE_CLASS, mxREAL);
     data2 = (float *) mxGetData(plhs[0]);
    
     // step
     plhs[1] = mxCreateNumericArray(1, dim2, mxINT32_CLASS, mxREAL);
     step = (int *) mxGetData(plhs[1]);

     // time
     plhs[2] = mxCreateNumericArray(1, dim2, mxSINGLE_CLASS, mxREAL);
     time = (float *) mxGetData(plhs[2]);

     // precision
     plhs[3] = mxCreateNumericArray(1, dim2, mxSINGLE_CLASS, mxREAL);
     precision = (float *) mxGetData(plhs[3]);
 
     // box    
     plhs[4] = mxCreateNumericArray(2, dim3, mxSINGLE_CLASS, mxREAL);
     box = (float *) mxGetData(plhs[4]);
   }
       
	if(mf->f == NULL)
		printf("File could not be opened! \n");
	do {
		breakVal = xtc_timestep(mf, ts, frame, data2, counter, numberAtoms, step, time, precision, box, startFrame, showmaxframes);
		if(breakVal == 1995) break;
		frame++;
	}while(frame < (showmaxframes));	   
}


