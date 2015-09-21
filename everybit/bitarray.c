/**
 * Copyright (c) 2012 MIT License by 6.172 Staff
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 **/

// Implements the ADT specified in bitarray.h as a packed array of bits; a bit
// array containing bit_sz bits will consume roughly bit_sz/8 bytes of
// memory.


#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include <sys/types.h>

#include "./bitarray.h"
#include "./byte_reverse.h"


//#define DEBUG

#ifdef DEBUG
  #define PRINT(VAR) printf("    "#VAR": %d\n", (int)(VAR));
#else
  #define PRINT(VAR)
#endif


// ********************************* Types **********************************

// Concrete data type representing an array of bits.
struct bitarray {
  // The number of bits represented by this bit array.
  // Need not be divisible by 8.
  size_t bit_sz;

  // The underlying memory buffer that stores the bits in
  // packed form (8 per byte).
  unsigned char *buf;
};

// Portable modulo operation that supports negative dividends.
//
// Many programming languages define modulo in a manner incompatible with its
// widely-accepted mathematical definition.
// http://stackoverflow.com/questions/1907565/c-python-different-behaviour-of-the-modulo-operation
// provides details; in particular, C's modulo
// operator (which the standard calls a "remainder" operator) yields a result
// signed identically to the dividend e.g., -1 % 10 yields -1.
// This is obviously unacceptable for a function which returns size_t, so we
// define our own.
//
// n is the dividend and m is the divisor
//
// Returns a positive integer r = n (mod m), in the range
// 0 <= r < m.
static size_t modulo(const ssize_t n, const size_t m);


// ******************************* Functions ********************************

bitarray_t *bitarray_new(const size_t bit_sz) {
  // Allocate an underlying buffer of ceil(bit_sz/8) bytes.
  unsigned char *const buf = calloc(1, bit_sz / 8 + ((bit_sz % 8 == 0) ? 0 : 1));
  if (buf == NULL) {
    return NULL;
  }

  // Allocate space for the struct.
  bitarray_t *const bitarray = malloc(sizeof(struct bitarray));
  if (bitarray == NULL) {
    free(buf);
    return NULL;
  }

  bitarray->buf = buf;
  bitarray->bit_sz = bit_sz;
  return bitarray;
}

void bitarray_free(bitarray_t *const bitarray) {
  if (bitarray == NULL) {
    return;
  }
  free(bitarray->buf);
  bitarray->buf = NULL;
  free(bitarray);
}

size_t bitarray_get_bit_sz(const bitarray_t *const bitarray) {
  return bitarray->bit_sz;
}

bool bitarray_get(const bitarray_t *const bitarray,
                  const size_t bit_index) {
  assert(bit_index < bitarray->bit_sz);

  unsigned char pos  = 7 - bit_index % 8;
  unsigned char mask = 1 << pos;
  return bitarray->buf[bit_index / 8] & mask;
}

void bitarray_set(bitarray_t *const bitarray,
                  const size_t bit_index,
                  const bool value) {
  assert(bit_index < bitarray->bit_sz);

  unsigned char pos  = 7 - bit_index % 8;
  unsigned char mask = ~(1 << pos);
  bitarray->buf[bit_index / 8] =
    (bitarray->buf[bit_index / 8] & mask) | (value << pos);
}

static inline void bitarray_set_byterange(bitarray_t *const ba,
                                          size_t bit_off,
                                          const size_t bit_len,
                                          const unsigned char val) {
  assert(ba);
  assert(bit_len <= 8);
  assert(bit_off + bit_len < ba->bit_sz);

  // TODO: possibly increase performance by setting two partial bytes?
  unsigned char mask = 1 << (bit_len - 1);
  for(size_t i = 0; i < bit_len; i++, bit_off++, mask >>= 1) {
     bitarray_set(ba, bit_off, val & mask);
  }
}

static inline void bitarray_reverse_byte(bitarray_t *const ba,
                                         size_t bit_off,
                                         size_t bit_len) {
  assert(ba);
  assert(bit_off + bit_len <= ba->bit_sz);  // ensure valid substring
  assert(bit_off % 8 + bit_len <= 8);       // ensure single byte

  size_t l = bit_off;
  size_t r = bit_off + bit_len - 1;
  bool tmp;
  
  // TODO: possibly pass in byte pointer, since it doesn't change?
  while(l < r) {
    tmp = bitarray_get(ba, l);
    bitarray_set(ba, l, bitarray_get(ba, r));
    bitarray_set(ba, r, tmp);
    l++;
    r--;
  }
}

static inline void bitarray_reverse_bytes(bitarray_t *const ba,
                                          unsigned char *left,
                                          unsigned char *right) {
  //TODO: possibly reverse all bits first, then swap bytes?
  while(left < right) {
    byte_reverse(left);
    byte_reverse(right);
    byte_swap(left, right);
    left++;
    right--;
  }

  // Reverse middle byte if odd number of bytes
  if(left == right) {
    PRINT(*left)
    byte_reverse(left);
    PRINT(*left)
  }
}

static inline void bitarray_shift_bytes(bitarray_t *const ba,
                                        unsigned char *left,
                                        unsigned char *right,
                                        signed char shift) {
  assert(ba);
  assert(left < right);

  unsigned char carry_mask;
  unsigned char carry_shift;
  unsigned char carry = 0;
  unsigned char tmp;

  if(shift > 0) {
    carry_shift = 8 - shift;
    carry_mask = (1 << shift) - 1;

    // Loop from left to right, shifting bits right in each byte
    for(; left <= right; left++) {
      tmp = (*left & carry_mask) << carry_shift;
      *left = carry | (*left >> shift);
      carry = tmp;
    }

    // Push carry bits into right edge byte
    *left = carry | (*left & (255 >> shift));
  } else if(shift < 0) {
    shift = -shift;
    carry_shift = 8 - shift;
    carry_mask = 255 << carry_shift;
    PRINT(shift)
    PRINT(carry_shift)
    PRINT(carry_mask)

    // Loop from right to left, shifting bits left in each byte
    for(; left <= right; right--) {
      tmp = (*right & carry_mask) >> carry_shift;
      *right = (*right << shift) | carry;
      carry = tmp;
      PRINT(*right)
      PRINT(carry)
    }

    // Push carry bits int left edge byte
    *right = (*right & (255 << shift)) | carry;
  }
}

static inline void bitarray_reverse(bitarray_t *const ba,
                                    size_t bit_off,
                                    size_t bit_len) {
  assert(ba);
  assert(bit_off + bit_len <= ba->bit_sz);  // ensure valid substring

  // Shortcut for single-byte range
  if(bit_off % 8 + bit_len <= 8) {
    // TODO: shortcut for complete byte using byte_reverse()?
    return bitarray_reverse_byte(ba, bit_off, bit_len);
  }

  printf("\n");
  PRINT(bit_off)
  PRINT(bit_len)

  // Determine length of left and right edge bits
  size_t left_edge_len  = 8 - (bit_off % 8);
  size_t right_edge_len = (bit_off + bit_len - 1) % 8 + 1;
  PRINT(left_edge_len)
  PRINT(right_edge_len)

  // Get first and last bytes
  unsigned char *first_byte = ba->buf + bit_off / 8;
  unsigned char *last_byte  = ba->buf + (bit_off + bit_len - 1) / 8;
  PRINT(*first_byte)
  PRINT(*last_byte)

  // Reverse bytes and bits in each byte
  bitarray_reverse_bytes(ba, first_byte + 1, last_byte - 1);

  // Retrieve edge bytes
  unsigned char left_edge = *first_byte & ~(255 << left_edge_len);
  unsigned char right_edge = *last_byte & ~(255 >> right_edge_len);
  PRINT(left_edge)
  PRINT(right_edge)

  // Reverse edge bytes
  byte_reverse(&left_edge);
  byte_reverse(&right_edge);
  PRINT(left_edge)
  PRINT(right_edge)

  // Calculate number of bits to shift in each byte
  signed char shift = right_edge_len - left_edge_len;
  PRINT(shift)

  // Shift bits in inner bytes
  bitarray_shift_bytes(ba, first_byte + 1, last_byte - 1, shift);

  // Swap edge bits
  left_edge  = left_edge >> (8 - left_edge_len);
  bitarray_set_byterange(ba, bit_off, right_edge_len, right_edge);
  bitarray_set_byterange(ba, bit_off + bit_len - left_edge_len, left_edge_len, left_edge);
}

void bitarray_rotate(bitarray_t *const bitarray,
                     const size_t bit_offset,
                     const size_t bit_length,
                     const ssize_t bit_right_amount) {
  assert(bit_offset + bit_length <= bitarray->bit_sz);

  // Nothing to rotate?
  if (bit_length == 0) return;

  // Get left rotation amount
  size_t n = modulo(-bit_right_amount, bit_length);

  // No rotation?
  if(n == 0) return;
  
  // Do the rotation
  bitarray_reverse(bitarray, bit_offset, n);
  bitarray_reverse(bitarray, bit_offset + n, bit_length - n);
  bitarray_reverse(bitarray, bit_offset, bit_length);
}

static size_t modulo(const ssize_t n, const size_t m) {
  const ssize_t signed_m = (ssize_t)m;
  assert(signed_m > 0);
  const ssize_t result = ((n % signed_m) + signed_m) % signed_m;
  assert(result >= 0);
  return (size_t)result;
}

