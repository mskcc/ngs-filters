PYTHON_HOME=/opt/common/CentOS_6-dev/python/python-2.7.10/bin/
export PATH=$PATH:/opt/common/CentOS_6-dev/bin/current/

    
    # Create virtualenv and install necessary packages
    /opt/common/CentOS_6-dev/python/python-2.7.10/bin/nosetests --with-xunit
#    pylint -f parseable myapp/ | tee pylint.out
