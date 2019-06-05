#ifndef PACK_COMMON_H
#define PACK_COMMON_H

#include <stdio.h>
#include <stdint.h>

inline void dumpPacket(uint8_t *p, uint32_t size, int n = 0)
{
    printf("data: \n");
    for(int i=0; i<size; i++){
      printf("%02x ",*(p+i));
      if( n>0 && (i+1)%n==0 )  printf("\n");
    }
    printf("\n");
}

#endif
