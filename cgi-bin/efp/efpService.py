#!/usr/bin/python3
"""
Module with classes for parsing service check configurations and handling the service checks
Created on Jan 5, 2010

@author: Robert Breit
"""
import http.client
import re
import socket
import urllib.error
import urllib.parse
import urllib.request
from xml.sax.handler import ContentHandler

import lxml.sax
from lxml import etree


class Service:
    def __init__(self, name, service_type):
        self.name = name
        self.type = service_type
        self.link = ''
        self.external = 'false'
        self.connect = ''
        self.icon = ''
        self.result_pattern = ''
        self.pattern_type = ''

    def add_connect(self, url):
        self.connect = url

    def add_icon(self, filename):
        self.icon = filename

    def add_link(self, url):
        self.link = url

    def add_external(self, webservice):
        self.external = webservice

    def add_no_result_regex(self, pattern):
        self.result_pattern = pattern
        self.pattern_type = 'negative'

    def add_result_regex(self, pattern):
        self.result_pattern = pattern
        self.pattern_type = 'positive'

    def check_service(self, gene):
        # Get sample signals through webservice.
        link = self.connect[:]
        if link == '':
            return 'Yes'  # return Yes, if no url defined

        link = re.sub(r"GENE", gene, link)

        try:
            # timeout condition
            timeout = 10
            socket.setdefaulttimeout(timeout)
            page = urllib.request.urlopen(link)
            result = page.read()
            check = re.search(self.result_pattern, result.decode('utf-8'))
            if self.pattern_type == 'negative':
                if check is None:
                    return 'Yes'  # result doesn't match pattern for NoResult
                else:
                    return None  # result matches pattern for NoResult
            else:
                if check is None:
                    return None  # result doesn't match pattern for Result
                else:
                    return 'Yes'  # result matches pattern for Result
        except (urllib.error.HTTPError, urllib.error.URLError, http.client.HTTPException):
            return 'error'

    def get_external(self):
        return self.external

    def get_link(self, gene):
        link = self.link[:]
        link = re.sub(r"GENE", gene, link)
        return link


class ServiceHandler(ContentHandler):
    def __init__(self, info):
        super().__init__()
        self.info = info
        self.currentService = ''

    def startElementNS(self, service_dict, qname, attrs):
        uri, name = service_dict
        if name == 'service':
            self.currentService = Service(attrs.getValueByQName('name'), attrs.getValueByQName('type'))

        if name == 'connect':
            self.currentService.add_connect(attrs.getValueByQName('url'))

        if name == 'icon':
            self.currentService.add_icon(attrs.getValueByQName('filename'))

        if name == 'link':
            self.currentService.add_link(attrs.getValueByQName('url'))

        if name == 'noresult_regex':
            self.currentService.add_no_result_regex(attrs.getValueByQName('pattern'))

        if name == 'result_regex':
            self.currentService.add_result_regex(attrs.getValueByQName('pattern'))

        if name == 'external':
            self.currentService.add_external(attrs.getValueByQName('webservice'))

    def endElementNS(self, qname, name):
        if name == 'service':
            self.info.add_service(self.currentService)


class Info:
    def __init__(self):
        self.services = {}  # Dictionary of views

    def add_service(self, service):
        self.services[service.name] = service

    def get_service(self, name):
        return self.services[name]

    def get_services(self):
        all_services = []
        external_services = []
        for name in self.services:
            service = self.get_service(name)
            if service.get_external() == 'true':
                external_services.append(name)
            else:
                all_services.append(name)
        all_services.extend(external_services)
        return all_services

    def load(self, file):
        handler = ServiceHandler(self)

        # Parse the file
        try:
            # Parse the file
            tree = etree.parse(file)
            lxml.sax.saxify(tree, handler)
        except etree.XMLSyntaxError:
            return 'error'

        return None
