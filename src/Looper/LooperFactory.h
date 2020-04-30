#ifndef LOOPERFACTORY_H
#define LOOPERFACTORY_H

#include <memory>
#include "Looper.h"

class LooperFactory {
public:
  static std::shared_ptr<Looper> create(const std::string &input_file);
};

#endif // LOOPERFACTORY_H
