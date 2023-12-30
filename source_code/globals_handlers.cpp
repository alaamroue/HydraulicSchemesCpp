/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include <iostream>
#include <CL/cl.h>

#define QUITE_SMALL     1E-9
#define	VERY_SMALL      1E-10
#define DOMAIN_DIR_N	0
#define DOMAIN_DIR_E	1
#define DOMAIN_DIR_S	2
#define DOMAIN_DIR_W	3
#define	GRAVITY         9.81



int get_group_id(int index) {
	return 1;
}
int get_global_size(int index) {
	return 1;
}
int get_global_id(int index) {
	return 1;
}
int get_local_id(int index) {
	return 1;
}
int get_local_size(int index) {
	return 1;
}