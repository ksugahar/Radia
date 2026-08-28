/* SPDX-License-Identifier: BSD-2-Clause
 * Shared TeX command-to-Unicode table.
 *
 * This is intentionally independent of the retired MTEF/.eqn codecs. It is
 * used by the TeX parser and structural-editing palette only.
 */
#ifndef RADIA_EQUATION_TEX_SYMBOLS_H
#define RADIA_EQUATION_TEX_SYMBOLS_H

#ifdef __cplusplus
extern "C" {
#endif

int tex_command_to_unicode(const char *cmd);
int tex_command_count(void);
const char *tex_command_name(int index);
int tex_command_code_at(int index);

#ifdef __cplusplus
}
#endif

#endif
