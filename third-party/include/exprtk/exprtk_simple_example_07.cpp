/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Simple Example 7                                           *
 * Author: Arash Partow (1999-2014)                           *
 * URL: http://www.partow.net/programming/exprtk/index.html   *
 *                                                            *
 * Copyright notice:                                          *
 * Free use of the Mathematical Expression Toolkit Library is *
 * permitted under the guidelines and in accordance with the  *
 * most current version of the Common Public License.         *
 * http://www.opensource.org/licenses/cpl1.0.php              *
 *                                                            *
 **************************************************************
*/


#include <cstdio>
#include <string>
#include "exprtk.hpp"


template <typename T>
void logic()
{
   typedef exprtk::expression<T> expression_t;
   std::string expression_string = "not(A and B) or C";

   exprtk::symbol_table<T> symbol_table;
   symbol_table.create_variable("A");
   symbol_table.create_variable("B");
   symbol_table.create_variable("C");

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   exprtk::parser<T> parser;
   parser.compile(expression_string,expression);

   printf(" # | A | B | C | %s\n"
          "---+---+---+---+-%s\n",
          expression_string.c_str(),
          std::string(expression_string.size(),'-').c_str());

   for (int i = 0; i < 8; ++i)
   {
      symbol_table.get_variable("A")->ref() = T(i & 0x01 ? 1 : 0);
      symbol_table.get_variable("B")->ref() = T(i & 0x02 ? 1 : 0);
      symbol_table.get_variable("C")->ref() = T(i & 0x04 ? 1 : 0);

      int result = static_cast<int>(expression.value());

      printf(" %d | %d | %d | %d | %d \n",
             i,
             static_cast<int>(symbol_table.get_variable("A")->value()),
             static_cast<int>(symbol_table.get_variable("B")->value()),
             static_cast<int>(symbol_table.get_variable("C")->value()),
             result);
   }
}

int main()
{
   logic<double>();
   return 0;
}
