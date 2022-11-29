//return date and time as text
void c_time_now(char text[]);

void c_time_now(char text[]){	
	char dbuffer[9];
	char tbuffer[9];
	_strdate(dbuffer);
	_strtime(tbuffer);
	strcpy(text, dbuffer);
	strcat(text, " & ");
	strcat(text, tbuffer);
}



