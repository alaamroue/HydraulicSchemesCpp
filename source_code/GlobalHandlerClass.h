/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include <iostream>

#pragma once
class GlobalHandlerClass {
public:
	int globalintX;
	int globalintY;

	GlobalHandlerClass(int, int);
	int get_group_id(int index);
	int get_global_size(int index);
	int get_global_id(int index);
	int get_local_id(int index);
	int get_local_size(int index);
};