#pragma once
#if KNL
#define _mm512_mullo_epi64(a, b) _mm512_mullo_epi32(a, b)
#define _mm256_maskz_loadu_epi32(mask,addr) _mm512_castsi512_si256(_mm512_maskz_loadu_epi32(mask,addr))
#define _mm256_mask_compressstoreu_epi32(addr,mask,a) _mm512_mask_compressstoreu_epi32(addr,mask,_mm512_castsi256_si512(a))
#define _mm256_cmpeq_epu32_mask(a,b) _mm512_cmpeq_epu32_mask(_mm512_castsi256_si512(a),_mm512_castsi256_si512(b))
#define _mm256_cmpgt_epu32_mask(a,b) _mm512_cmpgt_epu32_mask(_mm512_castsi256_si512(a),_mm512_castsi256_si512(b))
#define _mm256_mmask_i32gather_epi32(src,k,idx,addr,scale) _mm256_i32gather_epi32(addr,idx,scale)

#endif
