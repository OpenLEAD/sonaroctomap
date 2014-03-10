#include <boost/test/unit_test.hpp>
#include <sonaroctomap/Dummy.hpp>

using namespace sonaroctomap;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    sonaroctomap::DummyClass dummy;
    dummy.welcome();
}
