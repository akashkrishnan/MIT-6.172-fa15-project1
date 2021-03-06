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
#include <inttypes.h>

#include <sys/types.h>

#include "./bitarray.h"
#include "./byte_reverse.h"


// #define DEBUG

#ifdef DEBUG
  #define PRINT(VAR) printf("    "#VAR": %d (%02x)\n", (int)(VAR), (unsigned int)(VAR));
  #define PRINTP(VAR) printf("    "#VAR": %p\n", (VAR));
#else
  #define PRINT(VAR)
  #define PRINTP(VAR)
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

  unsigned char pos  = bit_index % 8;
  unsigned char mask = 1 << pos;
  return bitarray->buf[bit_index / 8] & mask;
}

void bitarray_set(bitarray_t *const bitarray,
                  const size_t bit_index,
                  const bool value) {
  assert(bit_index < bitarray->bit_sz);

  unsigned char pos  = bit_index % 8;
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
  assert(bit_off + bit_len <= ba->bit_sz);

  // TODO(akashk16): possibly increase performance by setting two partial bytes?
  for (size_t i = 0; i < bit_len; i++, bit_off++) {
    bitarray_set(ba, bit_off, val & (1 << i));
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

  // TODO(akashk16): possibly pass in byte pointer, since it doesn't change?
  // TODO(akashk16): we could use byte_reverse and masks to get the reversed bits
  while (l < r) {
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
  while (left < right) {
    byte_reverse_swap(left, right);
    left++;
    right--;
  }

  // Reverse middle byte if odd number of bytes
  if (left == right) {
    byte_reverse(left);
  }
}

static inline void bitarray_shift_bytes(bitarray_t *const ba,
                                        unsigned char *left,
                                        unsigned char *right,
                                        signed char shift) {
  assert(ba);

  uint64_t *cursor64;
  uint64_t carry_mask64;
  uint64_t carry_shift64;
  uint64_t carry64 = 0;
  uint64_t tmp64;

  unsigned char *cursor;
  unsigned char carry_mask;
  unsigned char carry_shift;
  unsigned char carry = 0;
  unsigned char tmp;

  if (shift > 0) {
    cursor64 = (uint64_t *)left;
    carry_shift64 = 64 - shift;
    carry_mask64  = 0xFFFFFFFFFFFFFFFF << carry_shift64;
    carry_shift   = 8 - shift;
    carry_mask    = 0xFF << carry_shift;

    while (cursor64 < (uint64_t *)(right - 8)) {
      tmp64 = (*cursor64 & carry_mask64) >> carry_shift64;
      *cursor64 = (*cursor64 << shift) | carry64;
      carry64 = tmp64;
      cursor64++;
    }

    carry = (unsigned char)carry64;
    cursor = (unsigned char *)cursor64;

    // Loop from left to right, shifting bits right (LE left shift) in each byte
    while (cursor <= right) {
      tmp = (*cursor & carry_mask) >> carry_shift;
      *cursor = (*cursor << shift) | carry;
      carry = tmp;
      cursor++;
    }

    // Push carry bits into right edge byte
    *cursor = (*cursor & (0xFF << shift)) | carry;
  } else if (shift < 0) {
    shift = -shift;
    cursor64 = (uint64_t *)(right - 7);
    carry_shift64 = 64 - shift;
    carry_mask64  = 0xFFFFFFFFFFFFFFFF >> carry_shift64;
    carry_shift   = 8 - shift;
    carry_mask    = 0xFF >> carry_shift;

    while (cursor64 >= (uint64_t *)left) {
      tmp64 = (*cursor64 & carry_mask64) << carry_shift64;
      *cursor64 = (*cursor64 >> shift) | carry64;
      carry64 = tmp64;
      cursor64--;
    }

    carry = (unsigned char)(carry64 >> 56);
    cursor = (unsigned char *)cursor64 + 7;

    // Loop from right to left, shifting bits left (LE left shift) in each byte
    while (cursor >= left) {
      tmp = (*cursor & carry_mask) << carry_shift;
      *cursor = (*cursor >> shift) | carry;
      carry = tmp;
      cursor--;
    }

    // Push carry bits into left edge byte
    *cursor = (*cursor & (0xFF >> shift)) | carry;
  }
}

static inline void bitarray_reverse(bitarray_t *const ba,
                                    size_t bit_off,
                                    size_t bit_len) {
  assert(ba);
  assert(bit_off + bit_len <= ba->bit_sz);  // ensure valid substring

  // Shortcut for single-byte range
  if (bit_off % 8 + bit_len <= 8) {
    // TODO(akashk16): shortcut for complete byte using byte_reverse()?
    return bitarray_reverse_byte(ba, bit_off, bit_len);
  }

  // Determine length of left and right edge bits
  size_t left_edge_len  = 8 - (bit_off % 8);
  size_t right_edge_len = (bit_off + bit_len - 1) % 8 + 1;

  // Get first and last bytes
  unsigned char *first_byte = ba->buf + bit_off / 8;
  unsigned char *last_byte  = ba->buf + (bit_off + bit_len - 1) / 8;

  // Reverse bytes and bits in each byte
  bitarray_reverse_bytes(ba, first_byte + 1, last_byte - 1);

  // Retrieve edge bytes
  unsigned char left_edge = *first_byte & ~(0xFF >> left_edge_len);
  unsigned char right_edge = *last_byte & ~(0xFF << right_edge_len);

  // Reverse edge bytes
  byte_reverse(&left_edge);
  byte_reverse(&right_edge);

  // Calculate number of bits to shift in each byte
  signed char shift = right_edge_len - left_edge_len;

  // Shift bits in inner bytes
  bitarray_shift_bytes(ba, first_byte + 1, last_byte - 1, shift);

  // Swap edge bits
  right_edge >>= (8 - right_edge_len);
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
  if (n == 0) return;

  // TODO(akashk16): add shortcuts for basic cases before doing the reverses?

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

